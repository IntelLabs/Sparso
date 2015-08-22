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
    int  matrix_version; // The latest version of the matrix. With version, we do not have to have another field "constant_valued".
    bool constant_structured; // A matrix might be changed, and thus its version updated, but its structure might stay the same
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
            delete schedule;
        }
    }

    void UpdateMatrixVersion(int new_version) {
        matrix_version = new_version;
        
        // Free all info based on the old matrix
        if (schedule != NULL) {
            delete schedule;
        }
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
        schedule = NULL;
    }
    
    ~BackwardTriangularSolveKnob() {
        if (schedule != NULL) {
            delete schedule;
        }
    }

    void UpdateMatrixVersion(int new_version) {
        matrix_version = new_version;
        
        // Free all info based on the old matrix
        if (schedule != NULL) {
            delete schedule;
        }
    }

private:
    // Info private to the backward triangular solver.
    int            matrix_version; // version of the matrix when the private info is built
    LevelSchedule* schedule;

    friend void BackwardTriangularSolve(
        int numrows, int numcols, int* colptr, int* rowval, double* nzval,
        double *y, const double *b, void* fknob);
};

/**************************** Usage of knobs *****************************/
void* NewMatrixKnob()
{
    MatrixKnob* m = new MatrixKnob;
    m->matrix_version = MIN_VALID_VERSION;
    m->constant_structured = false;
    return (void*)m;
}

void IncrementMatrixVersion(void* mknob)
{
    MatrixKnob* m = (MatrixKnob*)mknob;
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
    delete f;
}

void* NewBackwardTriangularSolveKnob()
{
    BackwardTriangularSolveKnob* f = new BackwardTriangularSolveKnob;
    return (void*)f;
}

void DeleteBackwardTriangularSolveKnob(void* fknob)
{
    BackwardTriangularSolveKnob* f = (BackwardTriangularSolveKnob*) fknob;
    delete f;
}

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
    CSR *A = new CSR(numrows, numcols, colptr, rowval, nzval, 1);
    LevelSchedule * schedule;
    if (fknob == NULL) {
        schedule = new LevelSchedule;
        schedule->constructTaskGraph(*A);
    } else {
        ForwardTriangularSolveKnob* f = (ForwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0];
        assert(m != NULL);
    
        if ((f->matrix_version == m->matrix_version) &&
            (f->schedule != NULL)) {
            schedule = f->schedule;
        } else {
            // Either no schedule, or is out of date
            f->UpdateMatrixVersion(m->matrix_version);
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(*A);
            f->schedule = schedule;
        }
    }
        
    int* invPerm = schedule->threadContToOrigPerm;
    forwardSolve(*A, y, b, *schedule, invPerm);
}

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob)
{
    CSR *A = new CSR(numrows, numcols, colptr, rowval, nzval, 1);
    LevelSchedule * schedule;
    if (fknob == NULL) {
        schedule = new LevelSchedule;
        schedule->constructTaskGraph(*A);
    } else {
        BackwardTriangularSolveKnob* f = (BackwardTriangularSolveKnob*) fknob;
        assert(f->mknobs.size() == 1);
        MatrixKnob* m = f->mknobs[0];
        assert(m != NULL);
    
        if ((f->matrix_version == m->matrix_version) &&
            (f->schedule != NULL)) {
            schedule = f->schedule;
        } else {
            // Either no schedule, or is out of date
            f->UpdateMatrixVersion(m->matrix_version);
            schedule = new LevelSchedule;
            schedule->constructTaskGraph(*A);
            f->schedule = schedule;
        }
    }
        
    int* invPerm = schedule->threadContToOrigPerm;
    backwardSolve(*A, y, b, *schedule, invPerm);
}
