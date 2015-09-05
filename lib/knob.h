#ifndef KNOB_H
#define KNOB_H

#ifdef __cplusplus
extern "C" {
#endif

/**************** Below is the interface for Julia *********************/

struct MatrixKnob;
struct FunctionKnob;

/**
 * colptr/rowval/nzval: CSC representation of sparse matrix
 */
MatrixKnob* NewMatrixKnob(int numrows, int numcols, int *colptr, int *rowval, double *nzval,
    bool constant_valued, bool constant_structured, bool is_symmetric, 
    bool is_structure_symmetric, bool is_structure_only);
void  SetConstantValued(MatrixKnob* mknob);
void  SetConstantStructured(MatrixKnob* mknob);
void  SetValueSymmetric(MatrixKnob* mknob);
void  SetStructureSymmetric(MatrixKnob *mknob);
void  SetStructureOnly(MatrixKnob *mknob);
void  SetMatrix(MatrixKnob* mknob, void* A);
void* GetDssHandle(MatrixKnob* mknob);
void* GetMatrix(MatrixKnob* mknob);
void  DeleteMatrixKnob(MatrixKnob* mknob);

/**
 * Compiler lets the library know if simple derivatives of the given matrix is already
 * available to avoid constructing them again.
 */
typedef enum
{
    DERIVATIVE_TYPE_TRANSPOSE = 0,
      /**< Derivative is transpose of the original matrix.
           For value derivative, in Julia notation, orig == derivative'
           For structure derivative, spones(orig) == spones(derivative') */
    DERIVATIVE_TYPE_SYMMETRIC = 1,
      /**< Derivative is a symmetric version of original matrix.
           For value derivative, original + original' - diag(original) == derivative
           For structure derivative, spones(original) + spones(original') == spones(derivative) */
    DERIVATIVE_TYPE_LOWER_TRIANGULAR = 2,
      /**< Derivative is lower triangular of original matrix
           For value derivative, tril(original) == derivative
           For structure derivative, spones(tril(original)) == spones(derivative)
           if A is lower triangular derivative of B, then B is symmetric derivative of A */
    DERIVATIVE_TYPE_UPPER_TRIANGULAR = 3,
      /**< Derivative is upper triangular of original matrix
           For value derivative, triu(original) == derivative
           For structure derivative, spones(triu(original)) == spones(derivative) 
           if A is upper triangular derivative of B, then B is symmetric derivative of B */
    DERIVATIVE_TYPE_COUNT,
} DerivativeType;

MatrixKnob* GetDerivative(MatrixKnob* mknob, DerivativeType type);
void  SetDerivative(MatrixKnob* mknob, DerivativeType type, MatrixKnob *derivative);

/**
 * Structural derivatives are where we are only interested in structures.
 * For example, SetStructureDerivative(A, DERIVATIVE_TYPE_TRANSPOSE, B) means
 * spones(A) == spones(B')
 */
MatrixKnob* GetStructureDerivative(MatrixKnob* mknob, DerivativeType type);
void  SetStructureDerivative(MatrixKnob* mknob, DerivativeType type, MatrixKnob *derivative);

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
    
FunctionKnob* NewFunctionKnob();
void DeleteFunctionKnob(FunctionKnob* fknob);

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
