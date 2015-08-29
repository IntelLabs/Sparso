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

// Function Knobs: each library function can choose to have a knob in order to take advantage of the context info in it.
// Otherwise, without a knob, the function just behaviors as usual.
class ForwardTriangularSolveKnob : public FunctionKnob {
public:
    ForwardTriangularSolveKnob() {
        matrix_version = INVALID_VERSION; // No valid private info built yet
        schedule = NULL;
    }
    
    ~ForwardTriangularSolveKnob() {
    }

    void UpdateMatrixVersion(int new_version) {
        matrix_version = new_version;
    }

private:
    // Info private to the forward triangular solver.
    int            matrix_version; // version of the matrix when the private info is built
    LevelSchedule* schedule;

    friend void ForwardTriangularSolve(
        int numrows, int numcols, int* colptr, int* rowval, double* nzval,
        double *y, const double *b, void* fknob);
};

class BackwardTriangularSolveKnob : public FunctionKnob {
public:
    BackwardTriangularSolveKnob() {
        matrix_version = INVALID_VERSION; // No valid private info built yet
    }
    
    ~BackwardTriangularSolveKnob() {
    }

    void UpdateMatrixVersion(int new_version) {
        matrix_version = new_version;
    }

private:
    // Info private to the backward triangular solver.
    int            matrix_version; // version of the matrix when the private info is built

    friend void BackwardTriangularSolve(
        int numrows, int numcols, int* colptr, int* rowval, double* nzval,
        double *y, const double *b, void* fknob);
};

class ADBKnob : public FunctionKnob {
public:
    ADBKnob() {
        structure = NULL;
    }
    
    ~ADBKnob() {
        if (structure != NULL) {
            delete (CSR *)structure;
        }
    }

private:
    CSR_Handle* structure;

    friend void* ADBInspect(const void *A, const void *B, void* fknob);
};

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
    DeleteOptimizedRepresentation(m);
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
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
    CSR *A = NULL;
    LevelSchedule *schedule;
    if (fknob == NULL) {
        A = new CSR(numrows, numcols, colptr, rowval, nzval, 1);
        A->computeInverseDiag();

        schedule = new LevelSchedule;
        schedule->constructTaskGraph(*A);
    } else {
        ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0];
        assert(m != NULL);
    
        if ((m->schedule != NULL) && 
            (m->constant_valued || m->constant_structured ||
             f->matrix_version == m->matrix_version)) {
            // FIXME: when matrix value has been changed (but not structure),
            // we need to update the values of shadow copy
            A = m->A;
            if (!A->idiag) A->computeInverseDiag();

            schedule = m->schedule;
        } else {
            CreateOptimizedRepresentation(m, numrows, numcols, colptr, rowval, nzval);
            A = m->A;
            A->computeInverseDiag();

            // Either no schedule, or is out of date
            f->UpdateMatrixVersion(m->matrix_version);
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(*A);
            m->schedule = schedule;
        }
    }
        
    int* invPerm = schedule->threadContToOrigPerm;
    forwardSolve(*A, y, b, *schedule, invPerm);

    if (!fknob) delete A;
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
    CSR *A = NULL;
    LevelSchedule * schedule;
    if (fknob == NULL) {
        A = new CSR(numrows, numcols, colptr, rowval, nzval, 1);
        A->computeInverseDiag();

        schedule = new LevelSchedule;
        schedule->constructTaskGraph(*A);
    } else {
        BackwardTriangularSolveKnob* f = (BackwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0];
        assert(m != NULL);
    
        if ((f->matrix_version == m->matrix_version) &&
            (m->schedule != NULL)) {
            A = m->A;
            if (!A->idiag) A->computeInverseDiag();

            schedule = m->schedule;
        } else {
            CreateOptimizedRepresentation(m, numrows, numcols, colptr, rowval, nzval);
            A = m->A;
            A->computeInverseDiag();
            
            // Either no schedule, or is out of date
            f->UpdateMatrixVersion(m->matrix_version);
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(*A);
            m->schedule = schedule;
        }
    }
        
    int* invPerm = schedule->threadContToOrigPerm;
    backwardSolve(*A, y, b, *schedule, invPerm);

    if (!fknob) delete A;
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

void* ADBInspect(const void *A, const void *B, void* fknob)
{
    assert(fknob != NULL);
    ADBKnob* f = (ADBKnob*) fknob;
    if (f->structure == NULL)
        f->structure = CSR_ADBInspect((CSR_Handle *)A, (CSR_Handle *)B);
    return (void*)(f->structure);
}
