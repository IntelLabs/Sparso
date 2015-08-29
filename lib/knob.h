#ifndef KNOB_H
#define KNOB_H

#ifdef __cplusplus
extern "C" {
#endif

/**************** Below is the interface for Julia *********************/

void* NewMatrixKnob();
void  IncrementMatrixVersion(void* mknob);
void* GetStructureProxy(void* mknob);
void* GetDssHandle(void* mknob);
void* GetMatrix(void* mknob);
void  DeleteMatrixKnob(void* mknob);

void  AddMatrixKnob(void* fknob, void* mknob);
void* GetMatrixKnob(void* fknob, int i);

void* NewForwardTriangularSolveKnob();
void  DeleteForwardTriangularSolveKnob(void* fknob);

void* NewBackwardTriangularSolveKnob();
void  DeleteBackwardTriangularSolveKnob(void* fknob);

void ForwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob);

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob);
    
void* NewADBKnob();

void DeleteADBKnob(void* fknob);

void *ADBInspect(
    const void *A, const void *B, void* fknob);

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
