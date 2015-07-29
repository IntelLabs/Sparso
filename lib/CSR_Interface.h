#ifndef CSR_INTERFACE_H
#define CSR_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CSR_Handle CSR_Handle;

// const/destruct
CSR_Handle *CSR_Create(int numRows, int numCols, int *rowptr, int *colidx, double *values, int base);
void CSR_Destroy(CSR_Handle *A);

// accessors
int CSR_GetNumRows(CSR_Handle *A);
int CSR_GetNumCols(CSR_Handle *A);
int CSR_GetNumNonZeros(CSR_Handle *A);

int *CSR_GetRowPtr(CSR_Handle *A);
int *CSR_GetColIdx(CSR_Handle *A);
double *CSR_GetValues(CSR_Handle *A);

// changing base
void CSR_Make0BasedIndexing(CSR_Handle *A);
void CSR_Make1BasedIndexing(CSR_Handle *A);

// load matrix
// Julia should call the following two function in order.
// Between the two calls, the CSR array space must be allocated.
void load_matrix_market_step1 (char *file, int *sizes, bool force_symmetric = false, bool transpose = false);
void load_matrix_market_step2(
  char *file, int *rowptr, int *colidx, double *values, int *sizes, bool one_based_CSR);

// C can directly call this once
void load_matrix_market(
  char *file,
  int **rowptr, int **colidx, double **values,
  int *is_symmetric, int *m, int *n, int *nnz,
  bool one_based_CSR = false, bool force_symmetric = false);

// w = alpha*A*x + beta*y + gamma
void CSR_MultiplyWithVector(
  double *w,
  double alpha, const CSR_Handle *A, const double *x,
  double beta, const double *y,
  double gamma);

// C = A*diag(d)*B
CSR_Handle *CSR_ADBInspect(
  const CSR_Handle *A, const CSR_Handle *B);
void CSR_ADB(
  CSR_Handle *C,
  const CSR_Handle *A, const CSR_Handle *B,
  const double *d);

void CSR_GetRCMPermutation(const CSR_Handle *A, int *perm, int *inversePerm);
void CSR_GetRCMPermutationWithoutPseudoDiameterSourceSelection(const CSR_Handle *A, int *perm, int *inversePerm);
void CSR_GetBFSPermutation(const CSR_Handle *A, int *perm, int *inversePerm);

void CSR_FindConnectedComponents(
  const CSR_Handle *A,
  int *numOfComponents, int **compToRoot, int **compSizes, int **compSizePrefixSum,
  int **nodesSortedByComp);

void CSR_Permute(const CSR_Handle *A, CSR_Handle *out, const int *columnPerm, const int *rowInversePerm);

int CSR_GetBandwidth(CSR_Handle *A);

void CSR_ReorderMatrix(int numRows, int numCols, int *rowptr, int *colidx, double *values, int *i1, int *j1, double *v1, 
                 int *perm, int *inversePerm, bool getPermutation, bool oneBasedInput, bool oneBasedOutput);
                 
// vector routines
void reorderVector(double *v, double *tmp, const int *perm, int len);
void reorderVectorWithInversePerm(double *v, double *tmp, const int *inversePerm, int len);

// w = alpha*x + beta*w
void waxpby(int n, double *w, double alpha, const double *x, double beta, const double *y);
// return dot(x*y)
double dot(int n, const double *x, const double *y);
// w = x./y
void pointwiseDivide(int n, double *w, const double *x, const double *y);
// w = x.*y
void pointwiseMultiply(int n, double *w, const double *x, const double *y);
// w = alpha*x + beta
void waxpb(int n, double *w, double alpha, const double *x, double beta);

#ifdef __cplusplus
}
#endif

#endif // CSR_INTERFACE_H
