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
    bool is_structure_symmetric, bool is_structure_only, bool is_single_def);
void  DeleteMatrixKnob(MatrixKnob* mknob);
void SetConstantValued(MatrixKnob* mknob);
bool IsConstantValued(MatrixKnob* mknob);
void SetConstantStructured(MatrixKnob* mknob);
bool IsConstantStructured(MatrixKnob* mknob);
void SetValueSymmetric(MatrixKnob *mknob);
bool IsValueSymmetric(MatrixKnob *mknob);
void SetStructureSymmetric(MatrixKnob *mknob);
bool IsStructureSymmetric(MatrixKnob *mknob);
void SetStructureOnly(MatrixKnob *mknob);
bool IsStructureOnly(MatrixKnob *mknob);
void SetMatrix(MatrixKnob* mknob, void* A);
void* GetMatrix(MatrixKnob* mknob);
void SetDssHandle(MatrixKnob* mknob, void* dss_handle);
void* GetDssHandle(MatrixKnob* mknob);
void PropagateMatrixInfo(MatrixKnob* to_mknob, MatrixKnob* from_mknob);

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

/**
 * w = alpha*A*x + beta*y + gamma
 */
void SpMV(
    int m, int n,
    double *w,
    double alpha,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double *x,
    double beta,
    double *y,
    double gamma,
    FunctionKnob *fknob);

void ForwardTriangularSolve(
    int L_numrows, int L_numcols, int* L_colptr, int* L_rowval, double* L_nzval,
    double *y, double *b, FunctionKnob *fknob);

void BackwardTriangularSolve(
    int U_numrows, int U_numcols, int* U_colptr, int* U_rowval, double* U_nzval,
    double *y, double *b, FunctionKnob *fknob);

void *CholFact(
    int m, int n, int *colptr, int *rowval, double *nzval,
    FunctionKnob *fknob);

void CholFactInverseDivide(
    void *factor, double *y, double *b, FunctionKnob *fknob);

// C = A*D*B
// A is m*k matrix
// B is k*n matrix
// C is m*n matrix
void ADB(
    int m, int n, int k,
    int **C_colptr, int **C_rowval, double **C_nzval,
    int *A_colptr, int *A_rowval, double *A_nzval,
    int *B_colptr, int *B_rowval, double *B_nzval,
    const double *d,
    FunctionKnob *fknob);

void ILU(
    int m, int n,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double *LU_nzval,
    FunctionKnob *fknob);
    
void IChol(
    int m, int n,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double *L_nzval,
    FunctionKnob *fknob);

void SpSquareWithEps(
    int m, int n,
    int **C_colptr, int **C_rowval, double **C_nzval,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double eps,
    FunctionKnob *fknob);

void SpAdd(
    int m, int n,
    int **C_colptr, int **C_rowval, double **C_nzval,
    double alpha,
    int *A_colptr, int *A_rowval, double *A_nzval,
    double beta,
    int *B_colptr, int *B_rowval, double *B_nzval,
    FunctionKnob *fknob);

double Trace(int n, int *A_colptr, int *A_rowval, double *A_nzval);

FunctionKnob* NewFunctionKnob();
void DeleteFunctionKnob(FunctionKnob* fknob);

/**
 * Let the function associated with fknob decide what permutation/inverse
 * permutation vector should be, and reorders its inputs and outputs accordingly.
 * Other functions just respect the decision, makes no decision, nor does any
 * reordering.
 */
void SetReorderingDecisionMaker(FunctionKnob *fknob);
int *GetRowReorderingVector(FunctionKnob *fknob, int *len);
int *GetRowInverseReorderingVector(FunctionKnob *fknob, int *len);
int *GetColReorderingVector(FunctionKnob *fknob, int *len);
int *GetColInverseReorderingVector(FunctionKnob *fknob, int *len);

/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif // KNOB_H
