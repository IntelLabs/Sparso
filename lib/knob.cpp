#include <stdio.h>
#include "knob.h"
#include <vector>
#include <assert.h>
#include "CSR_Interface.h"
#include "triSolve.hpp"
#include "SpMP/mm_io.h"
#include "SpMP/Utils.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"
#include <mkl.h>

using namespace SpMP;


// Julia variables' scope is function-wise at AST level, even though in the source
// level, it may appears to have nested scopes. Our Julia compiler creates matrix 
// knobs for array variables at the entry of a function, and deleted all
// of them at each exit of the function.
// Also, the Julia compiler creates function knobs for library function call sites at
// the entry of the function, associates it with the related matrix knobs, and 
// delete all of them at each exit of the function.

/**************************** Definition of knobs *****************************/

#define INVALID_VERSION   -1
#define MIN_VALID_VERSION 0

// A matrix knob stores the information shared by all the functions refering to this matrix.
// Function-specific info like a level schedule should be put into a function knob.
struct MatrixKnob {
    bool         constant_valued;     // The matrix is a constant in value(and thus of course constant in structure).
    int          matrix_version;      // The latest version of the matrix. It may be incremented dynamically if the matrix is not constant-valued.
    bool         constant_structured; // A matrix might be changed, and thus its version updated, but its structure might stay the same
    bool         is_symmetric;
    bool         is_structure_symmetric; // is_symmetric implies is_structure_symmetric
    bool         is_structure_only;   // true if structure of matrix is only used

    MatrixKnob  *derivatives[DERIVATIVE_TYPE_COUNT];

    // These are copied from Julia CSC matrix
    const int numrows;
    const int numcols;
    const int *colptr;
    const int *rowval;
    const double *nzval;

    MatrixKnob(int numrows, int numcols, const int *colptr, const int *rowval, const double *nzval) :
        numrows(numrows), numcols(numcols), colptr(colptr), rowval(rowval), nzval(nzval)
    {
    }

    // auxiliary data structure
    LevelSchedule     *schedule;
    _MKL_DSS_HANDLE_t dss_handle; 

    // The following fields shouldn't be passed to Julia because it has no idea
    // of their format.
    CSR        * A; // CSC used in Julia is often not be the best performing format. We want to decouple from it by having a shadow copy of the Julia CSC matrix in more optimized representation. For now, we fix it as CSR in SpMP. Later, we need to make it something like union so that we can change the matrix representation depending on the context
};

// This is the base class of all function knobs
struct FunctionKnob {
    std::vector<MatrixKnob *> mknobs; // knobs for the matrix arguments of the function
};

// Function Knobs: each library function can choose to have a function knob in 
// order to take advantage of the context info in it.
// Otherwise, without a knob, the function just behaviors as usual.
// The function knobs here do not have private data to save, and thus they are
// simply inheriting from the base class. 
class ForwardTriangularSolveKnob : public FunctionKnob { };
class BackwardTriangularSolveKnob : public FunctionKnob { };
class ADBKnob : public FunctionKnob { };
class CholfactKnob : public FunctionKnob { };
class CholmodFactorInverseDivideKnob : public FunctionKnob { };

/**************************** Usage of knobs *****************************/
// TODO: pass parameters (constant_structured, etc.) to NewMatrixKnob 
MatrixKnob* NewMatrixKnob(int numrows, int numcols, const int *colptr, const int *rowval, const double *nzval)
{
    MatrixKnob* m = new MatrixKnob(numrows, numcols, colptr, rowval, nzval);

    m->constant_valued = false;
    m->matrix_version = MIN_VALID_VERSION;
    m->constant_structured = false;
    m->is_symmetric = false;
    m->is_structure_symmetric = false;
    m->is_structure_only = false;

    for (int i = 0; i < DERIVATIVE_TYPE_COUNT; i++) {
        m->derivatives[i] = NULL;
    }

    m->schedule = NULL;
    m->A = NULL;

    return m;
}

static bool CheckMatrixKnobConsistency(MatrixKnob *m)
{
    if (m->constant_valued) {
        assert(m->constant_structured);
        if (!m->constant_structured) {
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

void IncrementMatrixVersion(MatrixKnob* mknob)
{
    assert(!mknob->constant_valued);
    mknob->matrix_version++;
}

void SetConstantValued(MatrixKnob* mknob)
{
    mknob->constant_valued = true;
    mknob->constant_structured = true;
}

void SetConstantStructured(MatrixKnob* mknob)
{
    mknob->constant_structured = true;
}

void SetValueSymmetric(MatrixKnob *mknob)
{
    mknob->is_symmetric = true;
    mknob->is_structure_symmetric = true;
}

void SetStructureSymmetric(MatrixKnob *mknob)
{
    mknob->is_structure_symmetric = true;
}

void SetStructureOnly(MatrixKnob *mknob)
{
    mknob->is_structure_only = true;
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

static CSR *CreateTransposedCSRCopyWith0BasedIndexing(
    int numrows, int numcols, const int *colptr, const int *rowval, const double *nzval)
{
    // assume input is 1-based indexing
    int nnz = colptr[numcols] - 1;
    int *colptr_copy = MALLOC(int, numcols + 1);
    int *rowval_copy = MALLOC(int, nnz);
    double *nzval_copy = NULL;
    if (nzval) {
        nzval_copy = MALLOC(double, nnz);
    }
    
    // Copy Julia CSC to SpMP CSR, while converting into 0
#pragma omp parallel for
    for (int i = 0; i <= numcols; i++)
        colptr_copy[i] = colptr[i] - 1;

    if (!nzval) {
#pragma omp parallel for
        for (int i = 0; i < nnz; i++) {
            rowval_copy[i] = rowval[i] - 1;
        }
    }
    else {
#pragma omp parallel for
        for (int i = 0; i < nnz; i++) {
            rowval_copy[i] = rowval[i] - 1;
            nzval_copy[i] = nzval[i];
        }
    }

    return new CSR(numcols, numrows, colptr_copy, rowval_copy, nzval_copy);
}

static void CreateOptimizedRepresentation(
    MatrixKnob *m, int numrows, int numcols, const int *colptr, const int *rowval, const double *nzval)
{
    // if something is not constant, it's a waste of time to create an optimized representation
    assert(m->is_structure_only && m->constant_structured || m->constant_valued);

    if (m->A) return;

    // Simply copying CSC to CSR will create a transposed version of original matrix
    CSR *AT = CreateTransposedCSRCopyWith0BasedIndexing(numrows, numcols, colptr, rowval, nzval);

    // When compiler tells you matrix must be symmetric, we skip checking symmetry
    bool needTranspose = false;
    if (!m->derivatives[DERIVATIVE_TYPE_TRANSPOSE]) {
        if (m->is_structure_only && !m->is_structure_symmetric) {
            needTranspose = !AT->isSymmetric(false);
        }
        if (!m->is_structure_only && !m->is_symmetric) {
            needTranspose = !AT->isSymmetric(true);
        }
    }
    if (needTranspose) {
        m->A = AT->transpose();

        MatrixKnob *knobTranspose = NewMatrixKnob(AT->m, AT->n, NULL, NULL, AT->values);
        m->derivatives[DERIVATIVE_TYPE_TRANSPOSE] = knobTranspose;
        knobTranspose->A = AT;
        knobTranspose->constant_valued = m->constant_valued;
        knobTranspose->constant_structured = m->constant_structured;
        knobTranspose->is_symmetric = m->is_symmetric;
        knobTranspose->is_structure_symmetric = m->is_structure_symmetric;
        knobTranspose->is_structure_only = m->is_structure_only;
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

MatrixKnob* GetDerivative(MatrixKnob* mknob, DerivativeType type)
{
    assert(type >= 0 && type < DERIVATIVE_TYPE_COUNT);
    return mknob->derivatives[type];
}

void SetDerivative(MatrixKnob *mknob, DerivativeType type, MatrixKnob *derivative)
{
    assert(type >= 0 && type < DERIVATIVE_TYPE_COUNT);
    if (DERIVATIVE_TYPE_TRANSPOSE == type || DERIVATIVE_TYPE_SYMMETRIC == type) {
        // no need to set transpose or symmetric version of A if A is symmetric
        assert(!mknob->is_symmetric);
        assert(!mknob->is_structure_only || !mknob->is_structure_symmetric);
    }
    mknob->derivatives[type] = derivative;
    if (derivative->constant_valued) {
    }
}

void* GetMatrix(MatrixKnob* mknob)
{
    return mknob->A;
}

void* GetDssHandle(MatrixKnob* mknob)
{
    return mknob->dss_handle;
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

FunctionKnob* NewForwardTriangularSolveKnob()
{
    return new ForwardTriangularSolveKnob;
}

void DeleteForwardTriangularSolveKnob(FunctionKnob* fknob)
{
    //delete fknob; FIXME - segfaults
}

FunctionKnob* NewBackwardTriangularSolveKnob()
{
    return new BackwardTriangularSolveKnob;
}

void DeleteBackwardTriangularSolveKnob(FunctionKnob* fknob)
{
    //delete fknob; FIXME - segfaults
}

//
// Forward triangular solve Ly = b, with A as the structure proxy of L.
// 
// ISSUES: 
// 1. So far, A is treated not only as a structure proxy, but also as
// as value proxy. This is not what we expected. However, if I want to use
// A just as a structure proxy, and use L itself for the value, the call to 
// forwardSolve(*L, y, b, *schedule, invPerm) throws segmentation fault. The
// reason is: L is not in the expected form to SpMP. It seems that to SpMP,
// a lower triangular matrix is represented with (1) non-diagonals, and (2) 
// idiag. However, in the incoming paramters L_colptr, L_rowval and L_nzval, 
// they include diagonals as well.
// 2. CSC vs. CSR (i.e. A vs. A'). Another hidden assumption of this function
// is it is forward solve for CSR array A. However, A is in CSC originally. So
// calling this function is actually forward solve for CSR array A'. 
// Now that A is not symmetric in value, but just in structure. Before calling
// this function, we have to compose A as (L + spdiagm(1./diag(A) * triu(A))'.
// Then we can call this function, and use A reliably as a value proxy of L.
// That is why I added the bandage in context-test2.jl pcg_symgs_with_context_opt().
// Please think how to get rid of this overhead. I think we have to use A only
// as structure proxy of L and U, NOT as value proxy. But then we have to solve 
// issue 1 above.
static void TriangularSolve(
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    double *y, const double *b, FunctionKnob* fknob,
    void (*solveFunction)(CSR&, double *, const double *, const LevelSchedule&, const int *))
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
            if (m->schedule) schedule = m->schedule;
            else {
                schedule = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC]->schedule;
                m->schedule = schedule;
            }
        } else if (m->constant_structured) {
            if (m->constant_valued) {
                CreateOptimizedRepresentation(m, L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
                L = m->A;
            }

            MatrixKnob *symKnob = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC];
            if (symKnob && (symKnob->A || symKnob->colptr)) {
                if (!symKnob->A) {
                    CreateOptimizedRepresentation(symKnob, L_numrows, L_numcols, symKnob->colptr, symKnob->rowval, symKnob->nzval);
                }
                schedule = new LevelSchedule;
                schedule->constructTaskGraph(*symKnob->A);
                m->derivatives[DERIVATIVE_TYPE_SYMMETRIC]->schedule = schedule;
                m->schedule = schedule;
            } else if (m->A) {
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

                    MatrixKnob *knobTranspose = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC];
                    if (!knobTranspose) {
                        knobTranspose = NewMatrixKnob(m->A->m, m->A->n, NULL, NULL, NULL);
                        m->derivatives[DERIVATIVE_TYPE_SYMMETRIC] = knobTranspose;
                        knobTranspose->constant_valued = m->constant_valued;
                        knobTranspose->constant_structured = m->constant_structured;
                        knobTranspose->is_symmetric = m->is_symmetric;
                        knobTranspose->is_structure_symmetric = m->is_structure_symmetric;
                        knobTranspose->is_structure_only = true;
                    }
                    knobTranspose->A = new CSR(m->A->m, m->A->n, symRowPtr, symColIdx, NULL);
                    knobTranspose->schedule = schedule;

                    FREE(symDiagPtr);
                    FREE(symExtPtr);
                }

                m->schedule = schedule;
            }
        } // constant_structured
    } // fknob != NULL

    if (!schedule && !L) {
        CSR *LT = CreateTransposedCSRCopyWith0BasedIndexing(
            L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
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
                MatrixKnob *knobTranspose = m->derivatives[DERIVATIVE_TYPE_SYMMETRIC];
                if (!knobTranspose) {
                    knobTranspose = NewMatrixKnob(m->A->m, m->A->n, NULL, NULL, NULL);
                    m->derivatives[DERIVATIVE_TYPE_SYMMETRIC] = knobTranspose;
                    knobTranspose->constant_valued = m->constant_valued;
                    knobTranspose->constant_structured = m->constant_structured;
                    knobTranspose->is_symmetric = m->is_symmetric;
                    knobTranspose->is_structure_symmetric = m->is_structure_symmetric;
                    knobTranspose->is_structure_only = true;
                }
                knobTranspose->A = new CSR(m->A->m, m->A->n, symRowPtr, symColIdx, NULL);
                knobTranspose->schedule = schedule;
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
            L = new CSR(L_numrows, L_numcols, L_colptr, L_rowval, L_nzval, 1);
            needToDeleteL = true;
        } else {
            CSR *LT = CreateTransposedCSRCopyWith0BasedIndexing(
                L_numrows, L_numcols, L_colptr, L_rowval, L_nzval);
            L = LT->transpose();
            delete LT;
            needToDeleteL = true;
        }
    }

    assert(L != NULL);
    assert(schedule != NULL);
    
    int* invPerm = schedule->threadContToOrigPerm;
    (*solveFunction)(*L, y, b, *schedule, invPerm);

    if (needToDeleteSchedule) {
        delete schedule;
    }
    if (needToDeleteL) {
        delete L;
    }
}

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, FunctionKnob* fknob)
{
    TriangularSolve(
        numrows, numcols, colptr, rowval, nzval, y, b, fknob,
        &forwardSolve);
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, FunctionKnob* fknob)
{
    TriangularSolve(
        numrows, numcols, colptr, rowval, nzval, y, b, fknob,
        &backwardSolve);
}

FunctionKnob* NewADBKnob()
{
    return new ADBKnob;
}

void DeleteADBKnob(FunctionKnob* fknob)
{
    delete fknob;
}

FunctionKnob* NewCholfactKnob()
{
    return new CholfactKnob;
}

void DeleteCholfactKnob(FunctionKnob* fknob)
{
    delete fknob;
}

FunctionKnob* NewCholmodFactorInverseDivideKnob()
{
    return new CholmodFactorInverseDivideKnob;
}

void DeleteCholmodFactorInverseDivideKnob(void* fknob)
{
    delete fknob;
}

/* vim: set tabstop=8 softtabstop=4 sw=4 expandtab: */
