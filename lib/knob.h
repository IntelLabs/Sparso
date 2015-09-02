#ifndef KNOB_H
#define KNOB_H

#ifdef __cplusplus
extern "C" {
#endif

/**************** Below is the interface for Julia *********************/

struct MatrixKnob;
struct FunctionKnob;

MatrixKnob* NewMatrixKnob(int numrows, int numcols, const int *colptr, const int *rowval, const double *nzval,
    bool constant_valued, bool constant_structured, bool is_symmetric, 
    bool is_structure_symmetric, bool is_structure_only);
void  IncrementMatrixVersion(MatrixKnob* mknob);
void  SetConstantValued(MatrixKnob* mknob);
void  SetConstantStructured(MatrixKnob* mknob);
void  SetValueSymmetric(MatrixKnob* mknob);
void  SetStructureSymmetric(MatrixKnob *mknob);
void  SetStructureOnly(MatrixKnob *mknob);
void* GetDssHandle(MatrixKnob* mknob);
void* GetMatrix(MatrixKnob* mknob);
void  DeleteMatrixKnob(MatrixKnob* mknob);

/**
 * Compiler lets the library knows if simple derivatives of the given matrix is already
 * available to avoid constructing them again.
 */
typedef enum
{
    DERIVATIVE_TYPE_TRANSPOSE = 0,
    DERIVATIVE_TYPE_SYMMETRIC = 1,
      /**< Proxy has sparsity structure of symmetric version of original matrix
           In Julia notation, orig == tril(derivative) && issym(derivative) */
    DERIVATIVE_TYPE_LOWER_TRIANGULAR = 2,
      /**< Proxy has sparsity structure of lower triangular of original matrix
           In Julia notation, tril(orig) == derivative */
    DERIVATIVE_TYPE_UPPER_TRIANGULAR = 3,
      /**< Proxy has sparsity structure of upper triangular of original matrix
           In Julia notation, triu(orig) == derivative */
    DERIVATIVE_TYPE_COUNT,
} DerivativeType;

MatrixKnob* GetDerivative(MatrixKnob* mknob, DerivativeType type);
void  SetDerivative(MatrixKnob* mknob, DerivativeType type, MatrixKnob *derivative);

void  AddMatrixKnob(FunctionKnob* fknob, MatrixKnob* mknob);
MatrixKnob* GetMatrixKnob(FunctionKnob* fknob, int i);

FunctionKnob* NewForwardTriangularSolveKnob();
void  DeleteForwardTriangularSolveKnob(FunctionKnob* fknob);

FunctionKnob* NewBackwardTriangularSolveKnob();
void  DeleteBackwardTriangularSolveKnob(FunctionKnob* fknob);

void ForwardTriangularSolve(
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    double *y, const double *b, FunctionKnob* fknob);

void BackwardTriangularSolve(
    int U_numrows, int U_numcols, int* U_colptr, int* U_rowval, double* U_nzval,
    double *y, const double *b, FunctionKnob* fknob);
    
FunctionKnob* NewADBKnob();
void  DeleteADBKnob(FunctionKnob* fknob);

FunctionKnob* NewCholfactKnob();
void  DeleteCholfactKnob(FunctionKnob* fknob);

FunctionKnob* NewCholmodFactorInverseDivideKnob();
void  DeleteCholmodFactorInverseDivideKnob(FunctionKnob* fknob);

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
