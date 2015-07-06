#include <stdlib.h>
#include "knob.h"

MatrixKnob* NewMatrixKnob(int m, int n, int* colptr, int* rowval, double* nzval)
{
    MatrixKnob* M = new MatrixKnob();
    M->version = 0;
    // M->representations[0]: // Jongsoo: please add whatever here with the input data
    M->representations[1] = NULL;
    return M;
}

void IncrementVersion(MatrixKnob* M)
{
    M->version++;
    // Jongsoo: please add free of memory of the representations.
    M->representations[0] = NULL; 
}

InverseDivideKnob* NewInverseDivideKnob(MatrixKnob* M)
{
    InverseDivideKnob* f = new InverseDivideKnob();
    f->matrix.M = M;
    f->matrix.version = M->version;
    f->schedule = NULL;
    return f;
}

double* InverseDivide(int m, int n, int* colptr, int* rowval, double* nzval, double* r, InverseDivideKnob* fknob)
{
    // Jongsoo: please add whatever with the data using the knob info.
}
