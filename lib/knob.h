#ifndef KNOB_H
#define KNOB_H

#ifdef __cplusplus
extern "C" {
#endif

/**************** Below is the interface for Julia *********************/

void* NewMatrixKnob();
void  IncrementMatrixVersion(void* mknob);
void  DeleteMatrixKnob(void* mknob);

void  AddMatrixKnob(void* fknob, void* mknob);

void* NewForwardTriangularSolveKnob();
void  DeleteForwardTriangularSolveKnob(void* fknob);

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob);

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
