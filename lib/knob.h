#ifndef KNOB_H
#define KNOB_H

#ifdef __cplusplus
extern "C" {
#endif

/**************** Below is the interface for Julia *********************/

void* NewMatrixKnob();
void  IncrementMatrixVersion(void* mknob);
void  SetConstantStructured(void* mknob);
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
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    int A_numrows, int A_numcols, int* A_colptr, int* A_rowval, double* A_nzval,
    double *y, const double *b, void* fknob);

void BackwardTriangularSolve(
    int numrows, int numcols, int* colptr, int* rowval, double* nzval,
    double *y, const double *b, void* fknob);
    
void* NewADBKnob();
void  DeleteADBKnob(void* fknob);

void* NewCholfactKnob();
void  DeleteCholfactKnob(void* fknob);

void* NewCholmodFactorInverseDivideKnob();
void  DeleteCholmodFactorInverseDivideKnob(void* fknob);

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
