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
    CSR        * structure_proxy;     // The structure of another matrix might represent the structure of this matrix. 

    CSR        * A; // CSC used in Julia is often not be the best performing format. We want to decouple from it by having a shadow copy of the Julia CSC matrix in more optimized representation. For now, we fix it as CSR in SpMP. Later, we need to make it something like union so that we can change the matrix representation depending on the context
    CSR        * AT; // to save transpose of A

    bool         is_symmetric;

    // auxiliary data structure
    LevelSchedule     *schedule;
    _MKL_DSS_HANDLE_t dss_handle; 
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
void* NewMatrixKnob()
{
    MatrixKnob* m = new MatrixKnob;
    m->constant_valued = false;
    m->matrix_version = MIN_VALID_VERSION;
    m->constant_structured = false;
    m->structure_proxy = NULL;
    m->A = NULL;
    m->AT = NULL;
    m->is_symmetric = false;
    m->schedule = NULL;
    return (void*)m;
}

void IncrementMatrixVersion(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    assert(!m->constant_valued);
    m->matrix_version++;
}

void SetConstantStructured(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    m->constant_structured = true;
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
    if (m->AT) {
      delete m->AT;
      m->AT = NULL;
    }
    m->is_symmetric = false;
}

static void CreateOptimizedRepresentation(
    MatrixKnob *m, int numrows, int numcols, int *colptr, int *rowval, double *nzval)
{
    DeleteOptimizedRepresentation(m);

    int nnz = colptr[numcols] - 1;
    int *colptr_copy = MALLOC(int, numcols + 1);
    int *rowval_copy = MALLOC(int, nnz);
    double *nzval_copy = MALLOC(double, nnz);
    
    // Copy Julia CSC to SpMP CSR, while converting into 0
#pragma omp parallel for
    for (int i = 0; i <= numcols; i++)
        colptr_copy[i] = colptr[i] - 1;

#pragma omp parallel for
    for (int i = 0; i < nnz; i++) {
        rowval_copy[i] = rowval[i] - 1;
        nzval_copy[i] = nzval[i];
    }

    // Simply copying CSC to CSR will create a transposed version of original matrix
    CSR *AT = new CSR(numcols, numrows, colptr_copy, rowval_copy, nzval_copy);

    if (!AT->isSymmetric(true, true)) {
      m->is_symmetric = false;

      m->AT = AT->transpose();
      m->AT = AT;
    }
    else {
      m->is_symmetric = true;
      m->A = AT;
    }
}

void* GetStructureProxy(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    return (void*) (m->structure_proxy);
}

void* GetMatrix(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    return (void*) (m->A);
}

void* GetDssHandle(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    return (void*) (m->dss_handle);
}

void DeleteMatrixKnob(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    //DeleteOptimizedRepresentation(m);
    if (m->schedule) {
        delete m->schedule;
        m->schedule = NULL;
    }
    delete m;
}

void AddMatrixKnob(void* fknob, void* mknob)
{
    assert(fknob != NULL);
    assert(mknob != NULL);
    FunctionKnob* f = (FunctionKnob*)fknob;
    MatrixKnob* m = (MatrixKnob*)mknob;
    f->mknobs.push_back(m);
}

void* GetMatrixKnob(void* fknob, int i)
{
    FunctionKnob* f = (FunctionKnob*)fknob;
    return (void*)(f->mknobs[i]);
}

void* NewForwardTriangularSolveKnob()
{
    ForwardTriangularSolveKnob* f = new ForwardTriangularSolveKnob;
    return (void*)f;
}

void DeleteForwardTriangularSolveKnob(void* fknob)
{
    ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
    //delete f; FIXME
}

void* NewBackwardTriangularSolveKnob()
{
    BackwardTriangularSolveKnob* f = new BackwardTriangularSolveKnob;
    return (void*)f;
}

void DeleteBackwardTriangularSolveKnob(void* fknob)
{
    BackwardTriangularSolveKnob* f = (BackwardTriangularSolveKnob*) fknob;
    //delete f; FIXME
}

void ForwardTriangularSolve(
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    int A_numrows, int A_numcols, int* A_colptr, int* A_rowval, double* A_nzval,
    double *y, const double *b, void* fknob)
{
/*    CSR * L = new CSR(L_numrows, L_numcols, L_colptr, L_rowval, L_nzval, 1);
  L->make0BasedIndexing();
    L->computeInverseDiag();
        for (int i = 0; i < L_numrows; i++) {
            for (int j = L->rowptr[i]; j < L->rowptr[i + 1]; j++) {
                if (i == L->colidx[j]) {
                    L->values[j] = 0.0;
                }
            }
        }
printf("L =");
L->print();
fflush(stdout);
    L->make1BasedIndexing();
*/
CSR* A = NULL;
    LevelSchedule* schedule = NULL;
    if (fknob == NULL) {
        A   = new CSR(A_numrows, A_numcols, A_colptr, A_rowval, A_nzval, 1);
        A->make0BasedIndexing();
        A->computeInverseDiag();
        schedule = new LevelSchedule;
        schedule->constructTaskGraph(*A);
        A->make1BasedIndexing();
    } else {
        ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0]; // This is the matrix knob for A
        assert(m != NULL);

        if ((m->schedule != NULL) && (m->constant_valued || m->constant_structured)) {
            schedule = m->schedule;
            A        = m->A;
        } else {
//        printf("reaching here\n");
        
            A   = new CSR(A_numrows, A_numcols, A_colptr, A_rowval, A_nzval, 1);
            A->make0BasedIndexing();
            A->computeInverseDiag();
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(*A);
            m->schedule = schedule;
//printf("A =");
//A->print();
//fflush(stdout);
            A->make1BasedIndexing();
            m->A = A;

//            int* invPerm = schedule->threadContToOrigPerm;
//    forwardSolve(*A, y, b, *schedule, invPerm);
//return;
        }
    }
/*
    assert(L != NULL);
    assert(schedule != NULL);
*/
    int* invPerm = schedule->threadContToOrigPerm;
    forwardSolve(*A, y, b, *schedule, invPerm);

    if (fknob == NULL) {
//        delete L;
        delete schedule;
    }
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
/*    CSR * L = new CSR(L_numrows, L_numcols, L_colptr, L_rowval, L_nzval, 1);
  L->make0BasedIndexing();
    L->computeInverseDiag();
        for (int i = 0; i < L_numrows; i++) {
            for (int j = L->rowptr[i]; j < L->rowptr[i + 1]; j++) {
                if (i == L->colidx[j]) {
                    L->values[j] = 0.0;
                }
            }
        }
printf("L =");
L->print();
fflush(stdout);
    L->make1BasedIndexing();
*/
    LevelSchedule* schedule = NULL;
    if (fknob == NULL) {
        CSR* A   = new CSR(numrows, numcols, colptr, rowval, nzval, 1);
        A->make0BasedIndexing();
        A->computeInverseDiag();
        schedule = new LevelSchedule;
        schedule->constructTaskGraph(*A);
        A->make1BasedIndexing();
    } else {
        ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0]; // This is the matrix knob for A
        assert(m != NULL);

        if ((m->schedule != NULL) && (m->constant_valued || m->constant_structured)) {
            schedule = m->schedule;
        } else {
        printf("reaching here\n");
        
            CSR* A   = new CSR(numrows, numcols, colptr, rowval, nzval, 1);
            A->make0BasedIndexing();
            A->computeInverseDiag();
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(*A);
            m->schedule = schedule;
//printf("A =");
//A->print();
//fflush(stdout);
            A->make1BasedIndexing();

            int* invPerm = schedule->threadContToOrigPerm;
    backwardSolve(*A, y, b, *schedule, invPerm);
return;
        }
    }
/*
    assert(L != NULL);
    assert(schedule != NULL);

    int* invPerm = schedule->threadContToOrigPerm;
    forwardSolve(*L, y, b, *schedule, invPerm);
*/
    if (fknob == NULL) {
//        delete L;
        delete schedule;
    }

}

void* NewADBKnob()
{
    ADBKnob* f = new ADBKnob;
    return (void*)f;
}

void DeleteADBKnob(void* fknob)
{
    ADBKnob* f = (ADBKnob*) fknob;
    delete f;
}

void* NewCholfactKnob()
{
    CholfactKnob* f = new CholfactKnob;
    return (void*)f;
}

void DeleteCholfactKnob(void* fknob)
{
    CholfactKnob* f = (CholfactKnob*) fknob;
    delete f;
}

void* NewCholmodFactorInverseDivideKnob()
{
    CholmodFactorInverseDivideKnob* f = new CholmodFactorInverseDivideKnob;
    return (void*)f;
}

void DeleteCholmodFactorInverseDivideKnob(void* fknob)
{
    CholmodFactorInverseDivideKnob* f = (CholmodFactorInverseDivideKnob*) fknob;
    delete f;
}


