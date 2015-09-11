#include <stdio.h>
#include "knob.h"
#include <vector>
#include <unordered_map>
#include <assert.h>
#include "CSR_Interface.h"
#include "TriSolve.hpp"
#include "SpMP/mm_io.h"
#include "SpMP/Utils.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"
#include <mkl.h>

using namespace std;
using namespace SpMP;

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
        perm(NULL), inverse_perm(NULL), numrows(numrows), numcols(numcols), colptr(colptr), rowval(rowval), nzval(nzval)
    {
    }

    // auxiliary data structure
    LevelSchedule     *schedule;
    _MKL_DSS_HANDLE_t dss_handle; 

    // The following fields shouldn't be passed to Julia because it has no idea
    // of their format.
    CSR        * A; // CSC used in Julia is often not be the best performing format. We want to decouple from it by having a shadow copy of the Julia CSC matrix in more optimized representation. For now, we fix it as CSR in SpMP. Later, we need to make it something like union so that we can change the matrix representation depending on the context

    int *perm, *inverse_perm;

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
    int *perm, *inverse_perm;
    int perm_len;

    FunctionKnob() : is_reordering_decision_maker(false), perm(NULL), inverse_perm(NULL) { }
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
    assert(constant_valued || constant_structured);

    assert(!constant_valued || constant_structured);
    assert(!is_symmetric || is_structure_symmetric);
    
    MatrixKnob* m = new MatrixKnob(numrows, numcols, colptr, rowval, nzval);

    m->constant_valued = constant_valued;
    m->constant_structured = constant_structured;
    m->is_symmetric = is_symmetric;
    m->is_structure_symmetric = is_structure_symmetric;
    m->is_structure_only = is_structure_only;
    m->is_single_def = is_single_def;

    for (int i = 0; i < DERIVATIVE_TYPE_COUNT; i++) {
        m->derivatives[i] = NULL;
        m->is_structure_derivatives[i] = false;
    }

    m->schedule = NULL;
    m->A = NULL;

    if (rowval) mknob_map[rowval] = m;

    return m;
}

void DeleteMatrixKnob(MatrixKnob* mknob)
{
    //DeleteOptimizedRepresentation(mknob);
    if (mknob->schedule) {
        delete mknob->schedule;
        mknob->schedule = NULL;
    }
    delete mknob;
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

    // Simply copying CSC to CSR will create a transposed version of original matrix
    CSR *AT = new CSR(numcols, numrows, colptr, rowval, nzval);

    // When compiler tells you matrix must be symmetric, we skip checking symmetry
    bool needTranspose = false;
    if (!m->derivatives[DERIVATIVE_TYPE_TRANSPOSE] || !m->is_structure_only && m->is_structure_derivatives[DERIVATIVE_TYPE_TRANSPOSE]) {
        if (m->is_structure_only && !m->is_structure_symmetric) {
            needTranspose = !AT->isSymmetric(false);
        }
        if (!m->is_structure_only && !m->is_symmetric) {
            needTranspose = !AT->isSymmetric(true);
        }
    }
    if (needTranspose) {
        m->A = AT->transpose();

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
    to_mknob->schedule   = from_mknob->schedule;
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
    delete fknob;
}

static void TriangularSolve(
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    double *y, double *b, FunctionKnob* fknob,
    void (*solveFunction)(CSR&, double *, const double *, const LevelSchedule&, const int *),
    void (*solveFunctionWithReorderedMatrix)(CSR&, double *, const double *, const LevelSchedule&))
{
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
            needToDeleteL = true;
        }
    }

    if (m && m->schedule && m->A) {
        fknob->perm = m->schedule->origToThreadContPerm;
        fknob->perm_len = m->A->m;
        fknob->inverse_perm = m->schedule->threadContToOrigPerm;

        if (m->perm != fknob->perm && fknob->is_reordering_decision_maker) {
            m->A = m->A->permute(fknob->perm, fknob->inverse_perm);
            m->perm = fknob->perm;
            m->inverse_perm = fknob->inverse_perm;

            reorderVectorWithInversePerm(b, fknob->inverse_perm, L->m);

            L = m->A;
        }
    }

    assert(L != NULL);
    assert(schedule != NULL);
    
    int* invPerm = schedule->threadContToOrigPerm;
    if (fknob->perm == m->perm) {
        (*solveFunctionWithReorderedMatrix)(*L, y, b, *schedule);
    }
    else {
        assert(m->perm == NULL);
        (*solveFunction)(*L, y, b, *schedule, invPerm);
    }

    if (needToDeleteSchedule) {
        delete schedule;
    }
    if (needToDeleteL) {
        delete L;
    }
}

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, double *b, FunctionKnob* fknob)
{
    TriangularSolve(
        numrows, numcols, colptr, rowval, nzval, y, b, fknob,
        &forwardSolve,
        &forwardSolveWithReorderedMatrix);
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, double *b, FunctionKnob* fknob)
{
    TriangularSolve(
        numrows, numcols, colptr, rowval, nzval, y, b, fknob,
        &backwardSolve,
        &backwardSolveWithReorderedMatrix);
}

void SetReorderingDecisionMaker(FunctionKnob *fknob)
{
    fknob->is_reordering_decision_maker = true;
}

int *GetReorderingVector(FunctionKnob *fknob, int *len)
{
    *len = fknob->perm_len;
    return fknob->perm;
}

int *GetInverseReorderingVector(FunctionKnob *fknob, int *len)
{
    *len = fknob->perm_len;
    return fknob->inverse_perm;
}

// The first few parameters (numRows to v) represent the source matrix to be reordered. The next few
// parameters (i1~v1) are the spaces that have been allocated to store the results. 
// Perm and inversePerm are the spaces that have been allocated for permutation and inverse permutation
// info; when getPermutation is true, this function computes and stores the info into them; otherwise,
// they already contain the info, and this function just uses it.
// oneBasedOutput: true if i1 and j1 should be 1-based indexing  
void CSR_ReorderMatrix(int numRows, int numCols, int *i, int *j, double *v, int *i1, int *j1, double *v1, 
                 int *perm, int *inversePerm, bool oneBasedOutput)
{
    double t1, t2, t3, t4, t5;
#ifdef PERF_TUNE
    t1 = omp_get_wtime();
#endif

    // The original and the result array space must be different
    assert(i != i1);
    assert(j != j1);    
    assert(v != v1);
    
    CSR *AT = new CSR(numRows, numCols, i, j, v);

#ifdef PERF_TUNE
    int orig_bw = AT->getBandwidth();
    t2 = omp_get_wtime();
#endif

#ifdef PERF_TUNE
    t3 = omp_get_wtime();
#endif

    CSR *newAT = new CSR(numRows, numCols, i1, j1, v1);
    AT->permuteRowptr(newAT, inversePerm);
    AT->permuteMain(newAT, perm, inversePerm);
    
#ifdef PERF_TUNE
    int rcm_bw = newAT->getBandwidth();
    t4 = omp_get_wtime();
#endif

    if (oneBasedOutput) {
        newAT->make1BasedIndexing();
    }
    else {
        newAT->make0BasedIndexing();
    }
    delete newAT;
    delete AT;

    auto itr = mknob_map.find(j);
    if (itr != mknob_map.end()) {
        itr->second->perm = perm;
        itr->second->inverse_perm = inversePerm;

        mknob_map[j1] = itr->second;

        itr->second->colptr = i1;
        itr->second->rowval = j1;
        itr->second->nzval = v1;
    }

#ifdef PERF_TUNE
        t5 = omp_get_wtime();
        double bytes = (double)i[numRows]*12;
        if (stats) {
          stats[0] += t5 - t1;
          stats[1] += t3 - t2;
          stats[2] += t4 - t3;
        }
    
#if 0
        printf("CSR_ReorderMatrix total: %f sec\n", t5 - t1);
        printf("\tCSR_GetRCMPermutation: %f sec (%f GB/s)\n", t3 - t2, bytes/(t3 - t2)/1e9);
        printf("\tCSR_Permute: %f sec (%f GB/s)\n", t4 - t3, bytes/(t4 - t3)/1e9);
        printf("\tBW changed: %d -> %d\n", orig_bw, rcm_bw);
        fflush(stdout);
#endif

#endif
}


/* vim: set tabstop=8 softtabstop=4 sw=4 expandtab: */
