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
    MatrixKnob * structure_proxy;     // The structure of another matrix might represent the structure of this matrix. 
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
        if (schedule != NULL) {
            //delete schedule;
        }
    }

    void UpdateMatrixVersion(int new_version) {
        matrix_version = new_version;
        
        // Free all info based on the old matrix
        if (schedule != NULL) {
            delete schedule;
            schedule = NULL;
        }
    }

private:
    // Info private to the forward triangular solver.
    int            matrix_version; // version of the matrix when the private info is built
    LevelSchedule* schedule;
    double*        idiag; // pre-computed inverse of diagonals

    friend void ForwardTriangularSolve(
        int numrows, int numcols, int* colptr, int* rowval, double* nzval,
        double *y, const double *b, void* fknob);
};

class BackwardTriangularSolveKnob : public FunctionKnob {
public:
    BackwardTriangularSolveKnob() {
        matrix_version = INVALID_VERSION; // No valid private info built yet
        schedule = NULL;
    }
    
    ~BackwardTriangularSolveKnob() {
        if (schedule != NULL) {
            //delete schedule;
        }
    }

    void UpdateMatrixVersion(int new_version) {
        matrix_version = new_version;
        
        // Free all info based on the old matrix
        if (schedule != NULL) {
            delete schedule;
            schedule = NULL;
        }
    }

private:
    // Info private to the backward triangular solver.
    int            matrix_version; // version of the matrix when the private info is built
    LevelSchedule* schedule;
    double*        idiag; // pre-computed inverse of diagonals

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
    return (void*)m;
}

void IncrementMatrixVersion(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
    assert(!m->constant_valued);
    m->matrix_version++;
}

void DeleteMatrixKnob(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
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

void* NewForwardTriangularSolveKnob()
{
    ForwardTriangularSolveKnob* f = new ForwardTriangularSolveKnob;
    return (void*)f;
}

void DeleteForwardTriangularSolveKnob(void* fknob)
{
    ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
    //delete f;
}

void* NewBackwardTriangularSolveKnob()
{
    BackwardTriangularSolveKnob* f = new BackwardTriangularSolveKnob;
    return (void*)f;
}

void DeleteBackwardTriangularSolveKnob(void* fknob)
{
    BackwardTriangularSolveKnob* f = (BackwardTriangularSolveKnob*) fknob;
    //delete f;
}

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
    CSR A(numrows, numcols, colptr, rowval, nzval, 1);
    LevelSchedule * schedule;
    if (fknob == NULL) {
        A.computeInverseDiag();

        schedule = new LevelSchedule;
        schedule->constructTaskGraph(A);
    } else {
        ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0];
        assert(m != NULL);
    
        if ((f->schedule != NULL) && 
            (m->constant_valued || m->constant_structured ||
             f->matrix_version == m->matrix_version)) {
            schedule = f->schedule;
            A.idiag = f->idiag;
        } else {
            A.computeInverseDiag();
            f->idiag = A.idiag;

            // Either no schedule, or is out of date
            f->UpdateMatrixVersion(m->matrix_version);
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(A);
            f->schedule = schedule;
        }
    }
        
    int* invPerm = schedule->threadContToOrigPerm;
    forwardSolve(A, y, b, *schedule, invPerm);

    A.idiag = NULL; // knob owns idiag
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
    CSR A(numrows, numcols, colptr, rowval, nzval, 1);
    LevelSchedule * schedule;
    if (fknob == NULL) {
        A.computeInverseDiag();

        schedule = new LevelSchedule;
        schedule->constructTaskGraph(A);
    } else {
        BackwardTriangularSolveKnob* f = (BackwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0];
        assert(m != NULL);
    
        if ((f->matrix_version == m->matrix_version) &&
            (f->schedule != NULL)) {
            schedule = f->schedule;
            A.idiag = f->idiag;
        } else {
            A.computeInverseDiag();
            f->idiag = A.idiag;
            
            // Either no schedule, or is out of date
            f->UpdateMatrixVersion(m->matrix_version);
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(A);
            f->schedule = schedule;
        }
    }
        
    int* invPerm = schedule->threadContToOrigPerm;
    backwardSolve(A, y, b, *schedule, invPerm);

    A.idiag = NULL; // knob owns idiag
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
