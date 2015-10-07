#include <cstdio>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <mkl.h>

#include "SpMP/mm_io.h"
#include "SpMP/Utils.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"

#include "knob.h"
#include "CSR_Interface.h"
#include "TriSolve.hpp"
#include "SpGEMM.hpp"
#include "ILU.hpp"
#include "IChol.hpp"
#include "BFSBipartite.h"

using namespace std;
using namespace SpMP;

struct ReorderingInfo {
    int *row_perm, *row_inverse_perm;
    int row_perm_len;
    int *col_perm, *col_inverse_perm;
    int col_perm_len;

    bool cost_benefit_analysis_done;

    ReorderingInfo() :
        row_perm(NULL), row_inverse_perm(NULL), row_perm_len(0),
        col_perm(NULL), col_inverse_perm(NULL), col_perm_len(0),
        cost_benefit_analysis_done(false)
    { }
};

static double bw_threshold_to_reorder = 45*1e9;
    // based on the result in endeavor, needs to be autotuned later

void SetBandwidthThresholdToReorder(double bw)
{
    bw_threshold_to_reorder = bw;
}

// set this true to see important decisions on matrix reordering
static bool LOG_REORDERING = false;
// set this true to see if we're doing transposes which is not cheap.
// we'd like to avoid transpose as much as possible if matrix
// is symmetric or a transposed version of the matrix is already
// available in the application context
static bool LOG_TRANSPOSE = false;

void SetLogLevel(int level)
{
    if (level > 0)
    {
        LOG_REORDERING = true;
        LOG_TRANSPOSE = true;
    }
    else
    {
        LOG_REORDERING = false;
        LOG_TRANSPOSE = false;
    }
}

static double spmp_spmv_time = 0;
static double spmp_trsv_time = 0;

static double knob_spmv_time = 0;
static double knob_trsv_time = 0;

double GetSpMPSpMVTime()
{
    return spmp_spmv_time;
}

void ResetSpMPSpMVTime()
{
    spmp_spmv_time = 0;
}

double GetSpMPTriangularSolveTime()
{
    return spmp_trsv_time;
}

void ResetSpMPTriangularSolveTime()
{
    spmp_trsv_time = 0;
}

double GetKnobSpMVTime()
{
    return knob_spmv_time;
}

void ResetKnobSpMVTime()
{
    knob_spmv_time = 0;
}

double GetKnobTriangularSolveTime()
{
    return knob_trsv_time;
}

void ResetKnobTriangularSolveTime()
{
    knob_trsv_time = 0;
}

// Julia variables' scope is function-wise at AST level, even though in the source
// level, it may appears to have nested scopes. Our Julia compiler creates matrix 
// knobs for array variables at the entry of a function, and deleted all
// of them at each exit of the function.
// Also, the Julia compiler creates function knobs for library function call sites at
// the entry of the function, associates it with the related matrix knobs, and 
// delete all of them at each exit of the function.

/**************************** Definition of knobs *****************************/

// A matrix knob stores the information shared by all the functions refering to this matrix.
// Function-specific info like a level schedule should be put into a function knob.
struct MatrixKnob {
    bool         constant_valued;     // The matrix is a constant in value(and thus of course constant in structure).
    bool         constant_structured; // A matrix might be changed, and thus its version updated, but its structure might stay the same
    bool         is_symmetric;
    bool         is_structure_symmetric; // is_symmetric implies is_structure_symmetric
    bool         is_structure_only;   // true if structure of matrix is only used
    bool         is_single_def; // The matrix is defined only once.

    MatrixKnob  *derivatives[DERIVATIVE_TYPE_COUNT];
    bool         is_structure_derivatives[DERIVATIVE_TYPE_COUNT];

    MatrixKnob(int numrows, int numcols, int *colptr, int *rowval, double *nzval) :
        schedule(NULL), dss_handle(NULL), A(NULL), numrows(numrows), numcols(numcols), colptr(colptr), rowval(rowval), nzval(nzval)
    {
        for (int i = 0; i < DERIVATIVE_TYPE_COUNT; i++) {
            derivatives[i] = NULL;
            is_structure_derivatives[i] = false;
        }
    }

    ~MatrixKnob() 
    {
        // Delete only the data owned by this matrix knob.
        // ISSUE: it seems that the fields in one matrix knob may be shared
        // by other matrix knobs. And deleting the fields causing segementation
        // fault when deleting the other matrix knobs.
        // TODO: maybe each knob should know which fields it owns, and which it
        // shares and does not own? Then each knob can safely releases the owned
        // fields. But we do not need this functionality for now.
        // if (schedule != NULL) delete schedule;
        // if (A != NULL) DeleteOptimizedRepresentation(this)
    }
    
    // auxiliary data structure
    LevelSchedule     *schedule;
    _MKL_DSS_HANDLE_t dss_handle; 

    // The following fields shouldn't be passed to Julia because it has no idea
    // of their format.
    CSR        * A; // CSC used in Julia is often not be the best performing format. We want to decouple from it by having a shadow copy of the Julia CSC matrix in more optimized representation. For now, we fix it as CSR in SpMP. Later, we need to make it something like union so that we can change the matrix representation depending on the context

    ReorderingInfo reordering_info;

    // These are copied from Julia CSC matrix
    int numrows;
    int numcols;
    int *colptr;
    int *rowval;
    double *nzval;
};

unordered_map<int *, MatrixKnob *> mknob_map; // rowval -> mknob map

// This is the base class of all function knobs
struct FunctionKnob {
    // knobs for the result and matrix arguments of the function, and maybe others
    // like consumers of the function. Each function knob may define this 
    // vector in its own way.
    std::vector<MatrixKnob *> mknobs;

    bool is_reordering_decision_maker;
    ReorderingInfo reordering_info;

    FunctionKnob() : is_reordering_decision_maker(false) { }
    ~FunctionKnob() 
    {
        // This does not affect the knobs themselves.
        mknobs.clear();
    } 
};

/**************************** Usage of knobs *****************************/
// TODO: pass parameters (constant_structured, etc.) to NewMatrixKnob 
MatrixKnob* NewMatrixKnob(int numrows, int numcols, int *colptr, int *rowval, double *nzval,
    bool constant_valued = false, bool constant_structured = false, bool is_symmetric = false, 
    bool is_structure_symmetric = false, bool is_structure_only = false,
    bool is_single_def = false)
{
    // The matrix knob is for a matrix that is constant either in value or
    // structure. Otherwise, it is not useful for now, although in future, we
    // might want to consider slowly changing matrices. 

    assert(!constant_valued || constant_structured);
    assert(!is_symmetric || is_structure_symmetric);
    
    MatrixKnob* m = new MatrixKnob(numrows, numcols, colptr, rowval, nzval);

    m->constant_valued = constant_valued;
    m->constant_structured = constant_structured;
    m->is_symmetric = is_symmetric;
    m->is_structure_symmetric = is_structure_symmetric;
    m->is_structure_only = is_structure_only;
    m->is_single_def = is_single_def;

    if (rowval) mknob_map[rowval] = m;

    return m;
}

void DeleteMatrixKnob(MatrixKnob* mknob)
{
    if (mknob != NULL)
    {
        mknob_map.erase(mknob->rowval);
        delete mknob;
    }
}

static bool CheckMatrixKnobConsistency(MatrixKnob *m)
{
    if (m->constant_valued) {
        assert(m->constant_structured);
        assert(!m->is_single_def); // There cannot be any definition.
        if (!m->constant_structured || m->is_single_def) {
            return false;
        }
    }
    if (m->is_symmetric) {
        assert(m->is_structure_symmetric);
        if (!m->is_structure_symmetric) {
            return false;
        }
    }

    return true;
}

void SetConstantValued(MatrixKnob* mknob)
{
    mknob->constant_valued = true;
    mknob->constant_structured = true;
}

bool IsConstantValued(MatrixKnob* mknob)
{
    return mknob->constant_valued;
}

void SetConstantStructured(MatrixKnob* mknob)
{
    mknob->constant_structured = true;
}

bool IsConstantStructured(MatrixKnob* mknob)
{
    return mknob->constant_structured;
}

void SetValueSymmetric(MatrixKnob *mknob)
{
    mknob->is_symmetric = true;
    mknob->is_structure_symmetric = true;
}

bool IsValueSymmetric(MatrixKnob *mknob)
{
    return mknob->is_symmetric;
}

void SetStructureSymmetric(MatrixKnob *mknob)
{
    mknob->is_structure_symmetric = true;
}

bool IsStructureSymmetric(MatrixKnob *mknob)
{
    return mknob->is_structure_symmetric;
}

void SetStructureOnly(MatrixKnob *mknob)
{
    mknob->is_structure_only = true;
}

bool IsStructureOnly(MatrixKnob *mknob)
{
    return mknob->is_structure_only;
}

void SetMatrix(MatrixKnob* mknob, void* A)
{
    mknob->A = (CSR*) A;
}

void* GetMatrix(MatrixKnob* mknob)
{
    return mknob->A;
}

void SetDssHandle(MatrixKnob* mknob, void* dss_handle)
{
    mknob->dss_handle = dss_handle;
}

void* GetDssHandle(MatrixKnob* mknob)
{
    return mknob->dss_handle;
}

static void DeleteOptimizedRepresentation(MatrixKnob *m)
{
    if (m->A) {
        FREE(m->A->rowptr);
        FREE(m->A->colidx);
        FREE(m->A->values);
        delete m->A;
        m->A = NULL;
    }
}

/**
 * nzval can be NULL, which means we're creating an optimized representation only for structure
 * of the matrix
 */
static void CreateOptimizedRepresentation(
    MatrixKnob *m, int numrows, int numcols, int *colptr, int *rowval, double *nzval)
{
    if (m->A) return;

    // if something is not constant, it's a waste of time to create an optimized representation
    if ((!m->is_structure_only || !m->constant_structured) && !m->constant_valued) {
        fprintf(stderr, "Warning: creating an optimized representation when not needed\n");
        return;
    }
    if (m->is_structure_only && nzval != NULL) {
        fprintf(stderr, "Warning: pass nzval as NULL when we only care about structure\n");
    }

    // If transpose is already available, just us it.
    MatrixKnob *knobTranspose = m->derivatives[DERIVATIVE_TYPE_TRANSPOSE];
    if (knobTranspose && (m->is_structure_only || !m->is_structure_derivatives[DERIVATIVE_TYPE_TRANSPOSE])) {
        m->A = new CSR(numrows, numcols, knobTranspose->colptr, knobTranspose->rowval, knobTranspose->nzval);
        return;
    }

    // Simply copying CSC to CSR will create a transposed version of original matrix
    CSR *AT = new CSR(numcols, numrows, colptr, rowval, nzval);

    // When compiler tells you matrix must be symmetric, we skip checking symmetry
    bool needTranspose = false;
    if (m->is_structure_only && !m->is_structure_symmetric) {
        needTranspose = !AT->isSymmetric(false);
    }
    if (!m->is_structure_only && !m->is_symmetric) {
        needTranspose = !AT->isSymmetric(true);
    }

    if (needTranspose) {
        m->A = AT->transpose();
        if (LOG_TRANSPOSE) {
            printf("Transposing a matrix\n");
        }

#ifndef NDEBUG
        AT->make0BasedIndexing();
        CSR *tempA = AT->transpose();
        tempA->make1BasedIndexing();
        AT->make1BasedIndexing();

        //m->A->print();
        //tempA->print();

        assert(m->A->equals(*tempA, true)); 
        delete tempA;
#endif

        MatrixKnob *knobTranspose = NewMatrixKnob(AT->m, AT->n, AT->rowptr, AT->colidx, AT->values,
                                                  m->constant_valued,
                                                  m->constant_structured,
                                                  m->is_symmetric,
                                                  m->is_structure_symmetric,
                                                  m->is_structure_only,
                                                  m->is_single_def);
        m->derivatives[DERIVATIVE_TYPE_TRANSPOSE] = knobTranspose;
        m->is_structure_derivatives[DERIVATIVE_TYPE_TRANSPOSE] = m->is_structure_only;
        knobTranspose->derivatives[DERIVATIVE_TYPE_TRANSPOSE] = m;
        knobTranspose->is_structure_derivatives[DERIVATIVE_TYPE_TRANSPOSE] = m->is_structure_only;
        knobTranspose->A = AT;
        assert(CheckMatrixKnobConsistency(knobTranspose));
    }
    else {
        m->A = AT;

        // in case the compiler was not able to prove the following facts
        m->is_structure_symmetric = true;
        if (!m->is_structure_only) {
            // in case the compiler was not able to prove the following facts
            m->is_symmetric = true;
        }
    }
}

void SetDerivative(MatrixKnob *mknob, DerivativeType type, MatrixKnob *derivative)
{
    assert(type >= 0 && type < DERIVATIVE_TYPE_COUNT);
    if (DERIVATIVE_TYPE_TRANSPOSE == type || DERIVATIVE_TYPE_SYMMETRIC == type) {
        // no need to set transpose or symmetric version of A if A is symmetric
        assert(!mknob->is_symmetric);
        assert(!mknob->is_structure_only || !mknob->is_structure_symmetric);
    }
    if (mknob->derivatives[type] && !mknob->is_structure_derivatives[type]) {
        fprintf(stderr, "Warning: mknob->derivative[%d] is already set.\n", type);
        return;
    }
    else {
        mknob->derivatives[type] = derivative;
        mknob->is_structure_derivatives[type] = false;
    }

    if (DERIVATIVE_TYPE_TRANSPOSE == type &&
            (!derivative->derivatives[type] || derivative->is_structure_derivatives[type])) {
        derivative->derivatives[type] = mknob;
        derivative->is_structure_derivatives[type] = false;
    }
    if ((DERIVATIVE_TYPE_LOWER_TRIANGULAR == type || DERIVATIVE_TYPE_UPPER_TRIANGULAR == type) &&
            (!derivative->derivatives[DERIVATIVE_TYPE_SYMMETRIC] || derivative->is_structure_derivatives[DERIVATIVE_TYPE_SYMMETRIC])) {
        derivative->derivatives[DERIVATIVE_TYPE_SYMMETRIC] = mknob;
        derivative->is_structure_derivatives[type] = false;
    }
}

MatrixKnob* GetDerivative(MatrixKnob* mknob, DerivativeType type)
{
    assert(type >= 0 && type < DERIVATIVE_TYPE_COUNT);
    return mknob->is_structure_derivatives[type] ? NULL : mknob->derivatives[type];
}

void SetStructureDerivative(MatrixKnob *mknob, DerivativeType type, MatrixKnob *derivative)
{
    assert(type >= 0 && type < DERIVATIVE_TYPE_COUNT);
    if (DERIVATIVE_TYPE_TRANSPOSE == type || DERIVATIVE_TYPE_SYMMETRIC == type) {
        // no need to set transpose or symmetric version of A if A is symmetric
        assert(!mknob->is_structure_symmetric);
    }
    if (mknob->derivatives[type]) {
        // we can overwrite structure derivative with value derivative,
        // but not vicex versa
        fprintf(stderr, "Warning: mknob->derivative[%d] is already set.\n", type);
        return;
    }
    else {
        mknob->derivatives[type] = derivative;
        mknob->is_structure_derivatives[type] = true;
    }

    if (DERIVATIVE_TYPE_TRANSPOSE == type && !derivative->derivatives[type]) {
        derivative->derivatives[type] = mknob;
        derivative->is_structure_derivatives[type] = true;
    }
    if ((DERIVATIVE_TYPE_LOWER_TRIANGULAR == type || DERIVATIVE_TYPE_UPPER_TRIANGULAR == type) && !derivative->derivatives[DERIVATIVE_TYPE_SYMMETRIC]) {
        derivative->derivatives[DERIVATIVE_TYPE_SYMMETRIC] = mknob;
        derivative->is_structure_derivatives[DERIVATIVE_TYPE_SYMMETRIC] = true;
    }
}

MatrixKnob *GetStructureDerivative(MatrixKnob *mknob, DerivativeType type)
{
    assert(type >= 0 && type < DERIVATIVE_TYPE_COUNT);
    return mknob->derivatives[type];    
}

void PropagateMatrixInfo(MatrixKnob* to_mknob, MatrixKnob* from_mknob)
{
    assert(to_mknob != NULL);
    assert(from_mknob != NULL);

    // Call to this function is due an assignment matrixA = matrixB. Thus
    // the destination matrix is not constant-valued. But both should be 
    // at least constant-structured, otherwise, they are useless.
    // The destination matrix can be a single-def, though, which means it 
    // is statically defined only by this source matrix, not anywhere else. In 
    // this situation, we can safely let the two matrices share all their
    // information. 
    assert(to_mknob->constant_structured && !to_mknob->constant_valued);
    assert(from_mknob->constant_structured || from_mknob->constant_valued);

    // Since the destination matrix is only constant-structured, we can only
    // copy the information that are determined only by structures, unless the
    // destination matrix is a single-def.

    // QUESTION: Schedule and dss_handle should depend only on structure, not value, right?
    if (from_mknob->schedule && !to_mknob->schedule)
        to_mknob->schedule   = from_mknob->schedule;
    if (from_mknob->dss_handle && !to_mknob->dss_handle)
        to_mknob->dss_handle = from_mknob->dss_handle;

    if (to_mknob->is_single_def) {
        to_mknob->A       = from_mknob->A;
        to_mknob->numrows = from_mknob->numrows;
        to_mknob->numcols = from_mknob->numcols;
        to_mknob->colptr  = from_mknob->colptr;
        to_mknob->rowval  = from_mknob->rowval;
        to_mknob->nzval   = from_mknob->nzval;
    }
    // How about derivatives?
    // Jongsoo: please see what to do here to be complete.
}

void AddMatrixKnob(FunctionKnob* fknob, MatrixKnob* mknob)
{
    assert(fknob != NULL);
    assert(mknob != NULL);
    fknob->mknobs.push_back(mknob);
}

MatrixKnob* GetMatrixKnob(FunctionKnob* fknob, int i)
{
    return fknob->mknobs[i];
}

FunctionKnob* NewFunctionKnob()
{
    return new FunctionKnob;
}

void DeleteFunctionKnob(FunctionKnob* fknob)
{
    if (fknob != NULL) delete fknob;
}

void SpMV(
    int m, int n,
    double *w,
    double alpha,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double *x,
    double beta,
    double *y,
    double gamma,
    FunctionKnob *fknob)
{
    knob_spmv_time -= omp_get_wtime();

    assert(fknob);
    MatrixKnob *mknob = GetMatrixKnob(fknob, 0);
    assert(mknob);

    CreateOptimizedRepresentation(mknob, m, n, A_colptr, A_rowval, A_nzval);

    double current_spmv_time = -omp_get_wtime();
    mknob->A->multiplyWithVector(w, alpha, x, beta, y, gamma);
    current_spmv_time += omp_get_wtime();
    spmp_spmv_time += current_spmv_time;

    if (fknob->is_reordering_decision_maker &&
            !fknob->reordering_info.cost_benefit_analysis_done &&
            fknob->reordering_info.row_inverse_perm == NULL) {

        fknob->reordering_info.cost_benefit_analysis_done = true;

        double bw = (12.*mknob->A->getNnz() + 8.*(mknob->A->m + mknob->A->n))/current_spmv_time;
        if (LOG_REORDERING) {
            printf("SpMV: reordering decision maker. BW measured is %g gbps\n", bw/1e9);
        }

        if (bw < bw_threshold_to_reorder) {
            double old_width_time = -omp_get_wtime();
            double old_w = mknob->A->getAverageWidth(true);
            old_width_time += omp_get_wtime();

            fknob->reordering_info.row_inverse_perm = MALLOC(int, m);
            fknob->reordering_info.col_perm = MALLOC(int, n);

            if (LOG_REORDERING) {
                printf("SpMV: permute matrix and input vector\n");
            }

            if (m != n) {
                CSR AT(n, m, A_colptr, A_rowval, A_nzval);

                fknob->reordering_info.row_perm = MALLOC(int, m);
                fknob->reordering_info.col_inverse_perm = MALLOC(int, n);

                double t = omp_get_wtime();
                bfsBipartite(
                    *mknob->A, AT,
                    fknob->reordering_info.row_perm,
                    fknob->reordering_info.row_inverse_perm,
                    fknob->reordering_info.col_perm,
                    fknob->reordering_info.col_inverse_perm);
                if (LOG_REORDERING) {
                    t = omp_get_wtime() - t;
                    printf("SpMV: matrix is rectangular. bfsBipartite takes %g (%g gbps)\n", t, ((double)mknob->A->getNnz()*12 + m*4*8)/t/1e9);
                }
            }
            else {
                double t = omp_get_wtime();
                mknob->A->getBFSPermutation(
                    fknob->reordering_info.col_perm,
                    fknob->reordering_info.row_inverse_perm);
                fknob->reordering_info.row_perm = fknob->reordering_info.col_perm;
                fknob->reordering_info.col_inverse_perm = fknob->reordering_info.row_inverse_perm;
                t = omp_get_wtime() - t;
                if (LOG_REORDERING) {
                    printf("SpMV: matrix is square. bfs takes %g (%g gbps)\n", t, ((double)mknob->A->getNnz()*12 + m*6*8)/t/1e9);
                }
            }

            double t = omp_get_wtime();

            CSR *newA = mknob->A->permute(
                fknob->reordering_info.col_perm,
                fknob->reordering_info.row_inverse_perm);

            double new_width_time = -omp_get_wtime();
            double new_w = newA->getAverageWidth();
            new_width_time += omp_get_wtime();

            if (LOG_REORDERING) {
                printf("SpMV: old_width = %g, new_width = %g, old_width takes %g sec new_width takes %g sec\n", old_w, new_w, old_width_time, new_width_time);
            }

            if (new_w/old_w < 1.2) {
                mknob->A = newA;

                mknob->reordering_info = fknob->reordering_info;

                reorderVectorWithInversePerm(
                    x, fknob->reordering_info.col_inverse_perm, n);
                if (y && y != x)
                    reorderVectorWithInversePerm(
                        y, fknob->reordering_info.row_perm, m);
                if (w != x)
                    reorderVectorWithInversePerm(
                        w, fknob->reordering_info.row_inverse_perm, n);

                if (LOG_REORDERING) {
                    printf("SpMV: row_perm=%p row_inverse_perm=%p col_perm=%p col_inverse_perm=%p\n", fknob->reordering_info.row_perm, fknob->reordering_info.row_inverse_perm, fknob->reordering_info.col_perm, fknob->reordering_info.col_inverse_perm);
                    t = omp_get_wtime() - t;
                    printf("SpMV: permutation takes %g (%g gbps)\n", t, ((double)mknob->A->getNnz()*2*12 + m*6*8)/t/1e9);
                }
            }
            else {
                if (fknob->reordering_info.row_perm != fknob->reordering_info.col_perm)
                    FREE(fknob->reordering_info.row_perm);
                fknob->reordering_info.row_perm = NULL;
                if (fknob->reordering_info.col_inverse_perm != fknob->reordering_info.row_inverse_perm)
                    FREE(fknob->reordering_info.col_inverse_perm);
                fknob->reordering_info.col_inverse_perm = NULL;
                FREE(fknob->reordering_info.col_perm);
                FREE(fknob->reordering_info.row_inverse_perm);

                delete newA;
            }
        }
    }

    knob_spmv_time += omp_get_wtime();
}

static void TriangularSolve_(
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    double *y, double *b, FunctionKnob* fknob,
    void (*solveFunction)(CSR&, double *, const double *, const LevelSchedule&),
    void (*solveFunctionWithReorderedMatrix)(CSR&, double *, const double *, const LevelSchedule&))
{
    knob_trsv_time -= omp_get_wtime();

    CSR* L = NULL;
    bool needToDeleteL = false;

    LevelSchedule* schedule = NULL;
    bool needToDeleteSchedule = false;

    MatrixKnob *m = NULL;

    if (fknob) {
        assert(fknob->mknobs.size() == 1);
        m = fknob->mknobs[0]; // This is the matrix knob for L or U
        assert(m != NULL);

        if ((m->schedule || m->derivatives[DERIVATIVE_TYPE_SYMMETRIC] && m->derivatives[DERIVATIVE_TYPE_SYMMETRIC]->schedule) && m->constant_structured) {
            // matrix has constant structure so reuse already available schedule
            if (m->schedule) schedule = m->schedule;
            else {
                schedule = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC]->schedule;
                m->schedule = schedule;
            }
        } else if (m->constant_structured) {
            // if schedule is not available construct it
            if (m->constant_valued) {
                // if it's constant valued, also save a shadow optimized representation
                CreateOptimizedRepresentation(m, L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
                L = m->A;
            }

            MatrixKnob *symKnob = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC];
            if (symKnob && (symKnob->A || symKnob->colptr)) {
                // if structurally symmetric version of this matrix is available use it for
                // constructing schedule
                CreateOptimizedRepresentation(symKnob, L_numrows, L_numcols, symKnob->colptr, symKnob->rowval, NULL);
                schedule = new LevelSchedule;
                schedule->constructTaskGraph(*symKnob->A);
                m->derivatives[DERIVATIVE_TYPE_SYMMETRIC]->schedule = schedule;
                m->schedule = schedule;
            } else if (m->A) {
                // if structurally symmetric version is not available, we have to construct it
                schedule = new LevelSchedule;

                int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
                bool wasSymmetric = getSymmetricNnzPattern(
                    m->A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);

                if (wasSymmetric) {
                    FREE(symRowPtr);
                    FREE(symColIdx);
                    FREE(symDiagPtr);
                    FREE(symExtPtr);

                    schedule->constructTaskGraph(*m->A);
                } else {
                    schedule->constructTaskGraph(
                       m->A->m,
                       symRowPtr, symDiagPtr, symExtPtr, symColIdx,
                       PrefixSumCostFunction(symRowPtr)); 

                    if (!symKnob) {
                        symKnob = NewMatrixKnob(m->A->m, m->A->n, NULL, NULL, NULL,
                                                m->constant_valued,
                                                m->constant_structured,
                                                m->is_symmetric,
                                                m->is_structure_symmetric,
                                                true /* is_structure_only */,
                                                m->is_single_def
                        );
                        m->derivatives[DERIVATIVE_TYPE_SYMMETRIC] = symKnob;
                        m->is_structure_derivatives[DERIVATIVE_TYPE_SYMMETRIC] = true;
                    }
                    symKnob->A = new CSR(m->A->m, m->A->n, symRowPtr, symColIdx, NULL);
                    symKnob->schedule = schedule;

                    FREE(symDiagPtr);
                    FREE(symExtPtr);
                }

                m->schedule = schedule;
            }
        } // constant_structured
    } // fknob != NULL

    if (!schedule && !L) {
        CSR *LT = new CSR(L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
        needToDeleteL = true;

        schedule = new LevelSchedule;
        needToDeleteSchedule = true;

        int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
        bool wasSymmetric = getSymmetricNnzPattern(
            LT, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);

        if (wasSymmetric) {
            FREE(symRowPtr);
            FREE(symColIdx);

            L = LT;
            schedule->constructTaskGraph(*L);

            if (m && m->constant_structured) {
                m->schedule = schedule;
                needToDeleteSchedule = false;
            }
        }
        else {
            schedule->constructTaskGraph(
               LT->m,
               symRowPtr, symDiagPtr, symExtPtr, symColIdx,
               PrefixSumCostFunction(symRowPtr)); 

            L = LT->transpose();
            delete LT;

            if (LOG_TRANSPOSE) {
                printf("Transposing a matrix\n");
            }

            if (m && m->constant_structured) {
                MatrixKnob *symKnob = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC];
                if (!symKnob) {
                    symKnob = NewMatrixKnob(m->A->m, m->A->n, NULL, NULL, NULL,
                                            m->constant_valued,
                                            m->constant_structured,
                                            m->is_symmetric,
                                            m->is_structure_symmetric,
                                            true /* is_structure_only */,
                                            m->is_single_def
                                            );
                    m->derivatives[DERIVATIVE_TYPE_SYMMETRIC] = symKnob;
                    m->is_structure_derivatives[DERIVATIVE_TYPE_SYMMETRIC] = true;
                }
                symKnob->A = new CSR(m->A->m, m->A->n, symRowPtr, symColIdx, NULL);
                symKnob->schedule = schedule;
                m->schedule = schedule;
                needToDeleteSchedule = false;
            } else {
                FREE(symRowPtr);
                FREE(symColIdx);
            }
        }

        FREE(symDiagPtr);
        FREE(symExtPtr);
    } else if (!L) {
        if (m->constant_valued) {
            CreateOptimizedRepresentation(m, L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
            L = m->A;
        } else if (m->is_symmetric) {
            L = new CSR(L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
            needToDeleteL = true;
        } else {
            CSR *LT = new CSR(L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
            L = LT->transpose();
            delete LT;

            if (LOG_TRANSPOSE) {
                printf("Transposing a matrix\n");
            }

            needToDeleteL = true;
        }
    }

    if (m && m->schedule && m->A) {
        fknob->reordering_info.col_perm = m->schedule->origToThreadContPerm;
        fknob->reordering_info.row_perm = m->schedule->origToThreadContPerm;
        fknob->reordering_info.row_inverse_perm = m->schedule->threadContToOrigPerm;
        fknob->reordering_info.col_inverse_perm = m->schedule->threadContToOrigPerm;
        fknob->reordering_info.row_perm_len = m->A->m;
        fknob->reordering_info.col_perm_len = m->A->n;
        assert(m->A->m == m->A->n);

        if (m->reordering_info.row_perm != fknob->reordering_info.row_perm &&
                !fknob->reordering_info.cost_benefit_analysis_done &&
                fknob->is_reordering_decision_maker) {

            if (LOG_REORDERING) {
                printf("TriangularSolve triSolve is reordering decision maker. Permute matrix and input vector\n");
            }
            if (LOG_REORDERING) {
                printf("TriangularSolve: row_perm=%p row_inverse_perm=%p col_perm=%p col_inverse_perm=%p\n", fknob->reordering_info.row_perm, fknob->reordering_info.row_inverse_perm, fknob->reordering_info.col_perm, fknob->reordering_info.col_inverse_perm);
            }

            fknob->reordering_info.cost_benefit_analysis_done = true;

            m->A = m->A->permute(
                fknob->reordering_info.col_perm,
                fknob->reordering_info.row_inverse_perm);

            m->reordering_info = fknob->reordering_info;

            reorderVectorWithInversePerm(
                b, fknob->reordering_info.col_inverse_perm, L->m);

            L = m->A;
        }
    }

    assert(L != NULL);
    assert(schedule != NULL);
    
    spmp_trsv_time -= omp_get_wtime();
    if (fknob->reordering_info.row_inverse_perm == m->reordering_info.row_inverse_perm) {
        (*solveFunctionWithReorderedMatrix)(*L, y, b, *schedule);
    }
    else {
        assert(m->reordering_info.row_inverse_perm == NULL);
        (*solveFunction)(*L, y, b, *schedule);
    }
    spmp_trsv_time += omp_get_wtime();

    if (needToDeleteSchedule) {
        delete schedule;
    }
    if (needToDeleteL) {
        delete L;
    }

    knob_spmv_time += omp_get_wtime();
}

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, double *b, FunctionKnob* fknob)
{
    TriangularSolve_(
        numrows, numcols, colptr, rowval, nzval, y, b, fknob,
        &forwardSolve,
        &forwardSolveWithReorderedMatrix);
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, double *b, FunctionKnob* fknob)
{
    TriangularSolve_(
        numrows, numcols, colptr, rowval, nzval, y, b, fknob,
        &backwardSolve,
        &backwardSolveWithReorderedMatrix);
}

static bool reuse_inspection = true;

void SetReuseInspection(bool reuse)
{
    reuse_inspection = reuse;
}

void *CholFact(
    int m, int n, int *colptr, int *rowval, double *nzval,
    FunctionKnob *fknob)
{
    assert(fknob);

    MatrixKnob *mknob_A = fknob->mknobs[0];
    assert(mknob_A);
    assert(mknob_A->constant_structured);

    if (!mknob_A->dss_handle) {
        MKL_INT opt = MKL_DSS_DEFAULTS;
        MKL_INT error = dss_create(mknob_A->dss_handle, opt);
        if (error != MKL_DSS_SUCCESS) {
            fprintf(stderr, "dss_create returned error code %d\n", error);
        }

        opt = MKL_DSS_NON_SYMMETRIC;
        int nnz = colptr[n] - 1;
        error = dss_define_structure(mknob_A->dss_handle, opt, colptr, m, m, rowval, nnz);
        if (error != MKL_DSS_SUCCESS) {
            fprintf(stderr, "dss_define_structure returned error code %d\n", error);
        }

        opt = MKL_DSS_AUTO_ORDER;
        error = dss_reorder(mknob_A->dss_handle, opt, NULL);
        if (error != MKL_DSS_SUCCESS) {
            fprintf(stderr, "dss_reorder returned error code %d\n", error);
        }
    }

    MKL_INT opt = MKL_DSS_DEFAULTS;
    MKL_INT error = dss_factor_real(mknob_A->dss_handle, opt, nzval);
    if (error != MKL_DSS_SUCCESS) {
        fprintf(stderr, "dss_factor_real returned error code %d\n", error);
    }

    void *ret= mknob_A->dss_handle;
    if (!reuse_inspection) {
        mknob_A->dss_handle = NULL;
    }
    return ret;
}

void CholFactInverseDivide(
    void *factor, double *y, double *b, FunctionKnob *fknob)
{
    assert(factor);
    _MKL_DSS_HANDLE_t dss_handle = (_MKL_DSS_HANDLE_t)factor;

    MKL_INT opt = MKL_DSS_DEFAULTS;
    int nrhs = 1;
    MKL_INT error = dss_solve_real(dss_handle, opt, b, nrhs, y);
    if (error != MKL_DSS_SUCCESS) {
        fprintf(stderr, "dss_solve_real returned error code %d\n", error);
    }
}

// C = A*D*B
// A is m*k matrix
// B is k*n matrix
// C is m*n matrix
void ADB(
    int m, int n, int k,
    int **C_colptr, int **C_rowval, double **C_nzval,
    int *A_colptr, int *A_rowval, double *A_nzval,
    int *B_colptr, int *B_rowval, double *B_nzval,
    const double *d,
    FunctionKnob *fknob)
{
    MatrixKnob *mknob_C = GetMatrixKnob(fknob, 0);
    MatrixKnob *mknob_A = GetMatrixKnob(fknob, 1);
    MatrixKnob *mknob_D = GetMatrixKnob(fknob, 2);
    MatrixKnob *mknob_B = GetMatrixKnob(fknob, 3);

    // By simply casting CSC as CSR, we get transposes of A and B.
    // Computing B'*D*A' gives (A*D*B)', so simply casting back the result
    // back as CSC will give the desired result, A*D*B.
    CSR *AT = new CSR(k, m, A_colptr, A_rowval, A_nzval);
    CSR *BT = new CSR(n, k, B_colptr, B_rowval, B_nzval);

    if (!mknob_C->A) {
        if (*C_colptr) {
            assert(*C_rowval);
            assert(*C_nzval);

            // Use memory passed from Julia, meaning
            // Julia uses in-place interface.
            mknob_C->A = new CSR(n, m, *C_colptr, *C_rowval, *C_nzval);
            inspectADB(mknob_C->A, BT, AT);
        }
        else {
            assert(!*C_rowval);
            assert(!*C_nzval);

            // Create memory inside
            mknob_C->A = inspectADB(BT, AT);
        }
    }

    adb(mknob_C->A, BT, AT, d);

    if (*C_colptr) {
        assert(*C_colptr == mknob_C->A->rowptr);
        assert(*C_rowval == mknob_C->A->colidx);
        assert(*C_nzval == mknob_C->A->values);
    }
    else {
        assert(!*C_rowval);
        assert(!*C_nzval);

        *C_colptr = mknob_C->A->rowptr;
        *C_rowval = mknob_C->A->colidx;
        *C_nzval = mknob_C->A->values;
    }

    if (!reuse_inspection) {
        mknob_C->A = NULL;
    }
}

void ILU(
    int m, int n,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double *LU_nzval,
    FunctionKnob *fknob)
{
    CSR *A = new CSR(m, n, A_colptr, A_rowval, A_nzval);
    LevelSchedule *schedule = new LevelSchedule;
    schedule->constructTaskGraph(*A);
    ilu0(LU_nzval, *A, *schedule);
    delete schedule;
    delete A;
}

void IChol(
    int m, int n,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double *L_nzval,
    FunctionKnob *fknob)
{
    CSR *A = new CSR(m, n, A_colptr, A_rowval, A_nzval);
    LevelSchedule *schedule = new LevelSchedule;
    schedule->constructTaskGraph(*A);
    ichol0(*A, L_nzval, *schedule);
    delete schedule;
    delete A;
}

void SpSquareWithEps(
    int m, int n,
    int **C_colptr, int **C_rowval, double **C_nzval,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double eps,
    FunctionKnob *fknob)
{
    CSR AT(n, m, A_colptr, A_rowval, A_nzval);

    CSR *A;
    if (fknob && fknob->mknobs[0] && fknob->mknobs[0]->is_symmetric) {
        A = &AT;
    }
    else {
        A = AT.transpose();

        if (LOG_TRANSPOSE) {
            printf("Transposing a matrix\n");
        }
    }

    CSR *C = SpGEMMWithEps(A, &AT, eps);

    if (A != &AT) {
        delete A;
    }

    *C_colptr = C->rowptr;
    *C_rowval = C->colidx;
    *C_nzval = C->values;
}

void SpAdd(
    int m, int n,
    int **C_colptr, int **C_rowval, double **C_nzval,
    double alpha,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double beta,
    int *B_colptr, int *B_rowval, double *B_nzval,
    FunctionKnob *fknob)
{
    CSR AT(n, m, A_colptr, A_rowval, A_nzval);
    CSR BT(n, m, B_colptr, B_rowval, B_nzval);
    // It doesn't matter if we do matrix addition in CSR or CSC
    CSR *C = SpAdd(alpha, &AT, beta, &BT);

    *C_colptr = C->rowptr;
    *C_rowval = C->colidx;
    *C_nzval = C->values;
}

double Trace(int n, int *A_colptr, int *A_rowval, double *A_nzval)
{
    double trace = 0;
#pragma omp parallel for reduction(+:trace)
    for (int i = 0; i < n; ++i) {
        for (int j = A_colptr[i] - 1; j < A_colptr[i + 1] - 1; ++j) {
            if (A_rowval[j] - 1 == i) {
                trace += A_nzval[j];
            }
        }
    }
    return trace;
}

void SetReorderingDecisionMaker(FunctionKnob *fknob)
{
    fknob->is_reordering_decision_maker = true;
}

int *GetRowReorderingVector(FunctionKnob *fknob, int *len)
{
    *len = fknob->reordering_info.row_perm_len;
    return fknob->reordering_info.row_perm;
}

int *GetRowInverseReorderingVector(FunctionKnob *fknob, int *len)
{
    *len = fknob->reordering_info.row_perm_len;
    return fknob->reordering_info.row_inverse_perm;
}

int *GetColReorderingVector(FunctionKnob *fknob, int *len)
{
    *len = fknob->reordering_info.col_perm_len;
    return fknob->reordering_info.col_perm;
}

int *GetColInverseReorderingVector(FunctionKnob *fknob, int *len)
{
    *len = fknob->reordering_info.col_perm_len;
    return fknob->reordering_info.col_inverse_perm;
}

void ReorderMatrixInplace(int numRows, int numCols, int *colptr, int *rowval, double *nzval, 
                 int *perm, int *inversePerm)
{
    assert(perm);
    assert(inversePerm);

    CSR AT(numCols, numRows, colptr, rowval, nzval);
    CSR newAT(AT);

    newAT.permuteRowptr(&AT, inversePerm);
    newAT.permuteMain(&AT, perm, inversePerm);
    
    if (LOG_REORDERING) {
        printf("CSR_ReorderMatrix: perm=%p inversePerm=%p\n", perm, inversePerm);
    }

    auto itr = mknob_map.find(rowval);
    if (itr != mknob_map.end()) {
        if (LOG_REORDERING) {
            printf("CSR_ReorderMatrix: a related matrix knob found. Set its reordering information\n");
        }

        itr->second->reordering_info.row_perm = perm;
        itr->second->reordering_info.col_perm = perm;
        itr->second->reordering_info.col_inverse_perm = inversePerm;
        itr->second->reordering_info.row_inverse_perm = inversePerm;
        itr->second->reordering_info.row_perm_len = numRows;
        itr->second->reordering_info.col_perm_len = numCols;

        if (itr->second->A) {
            if (LOG_REORDERING) {
                printf("CSR_ReorderMatrix: an optimized matrix representation found. Update it.\n");
            }

            itr->second->A = NULL;
            CreateOptimizedRepresentation(
                itr->second, numRows, numCols, colptr, rowval, nzval);
        }
    }
}

void ReorderMatrix(int numRows, int numCols, int *colptr, int *rowval, double *nzval, int *colptr_out, int *rowval_out, double *nzval_out, 
                 int *perm, int *inversePerm)
{
    // The original and the result array space must be different
    assert(colptr != colptr_out);
    assert(rowval != rowval_out);    
    assert(nzval != nzval_out);
    
    assert(perm);
    assert(inversePerm);
    
    CSR AT(numCols, numRows, colptr, rowval, nzval);
    CSR newAT(numCols, numRows, colptr_out, rowval_out, nzval_out);

    AT.permuteRowptr(&newAT, inversePerm);
    AT.permuteMain(&newAT, perm, inversePerm);

    if (LOG_REORDERING) {
        printf("CSR_ReorderMatrix: perm=%p inversePerm=%p\n", perm, inversePerm);
    }

    auto itr = mknob_map.find(rowval);
    if (itr != mknob_map.end()) {
        if (LOG_REORDERING) {
            printf("CSR_ReorderMatrix: a related matrix knob found. Set its reordering information\n");
        }

        itr->second->reordering_info.row_perm = perm;
        itr->second->reordering_info.col_perm = perm;
        itr->second->reordering_info.col_inverse_perm = inversePerm;
        itr->second->reordering_info.row_inverse_perm = inversePerm;
        itr->second->reordering_info.row_perm_len = numRows;
        itr->second->reordering_info.col_perm_len = numCols;

        mknob_map[rowval_out] = itr->second;

        itr->second->colptr = colptr_out;
        itr->second->rowval = rowval_out;
        itr->second->nzval = nzval_out;

        if (itr->second->A) {
            if (LOG_REORDERING) {
                printf("CSR_ReorderMatrix: an optimized matrix representation found. Update it.\n");
            }

            itr->second->A = NULL;
            CreateOptimizedRepresentation(
                itr->second, numRows, numCols, colptr_out, rowval_out, nzval_out);
        }
    }
}

/* vim: set tabstop=8 softtabstop=4 sw=4 expandtab: */
