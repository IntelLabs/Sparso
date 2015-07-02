#ifndef KNOB_H
#define KNOB_H

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_NUM_SPARSE_MATRIX_REPS 10

/*
typedef enum {
    SPARSE_MATRIX_CSC = 0,
    SPARSE_MATRIX_CSR = 1
} MatrixFormat;
    
struct MatrixRepresentation {
    bool         oneBased; // Is the representation based on 1? If not, it is based on 0
    MatrixFormat format;
    void*        representation; // Exactly how to represent and use the matrix is the job of the library
};
*/

struct MatrixKnob {
    int   version; // The latest version of the matrix.
    void* representations[MAX_NUM_SPARSE_MATRIX_REPS]; // Representations of the matrix with the latest version. Stop at NULL. Exactly how to represent and use the matrix is the job of the library
};

struct InverseDivideKnob {
    struct {
        MatrixKnob* M;
        int         version; // private version of the matrix
    } matrix;
    void* schedule; // Exactly how to represent and use the schedule (task graph included) is the job of the library 
};

/**************** Below is the interface for Julia *********************/

MatrixKnob*        NewMatrixKnob(int m, int n, int* colptr, int* rowval, double* nzval);
void               IncrementVersion(MatrixKnob* M); // Increment version. Clear all cached data with the old version.
InverseDivideKnob* NewInverseDivideKnob(MatrixKnob* M);
double*            InverseDivide(int m, int n, int* colptr, int* rowval, double* nzval, double* r, InverseDivideKnob* fknob = NULL); // M\r

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
