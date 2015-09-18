#ifndef CSR_INTERFACE_H
#define CSR_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CSR_Handle CSR_Handle;

// const/destruct
CSR_Handle *CSR_Create(int numRows, int numCols, int *rowptr, int *colidx, double *values);
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

void load_vector_matrix_market(const char *fileName, double **v, int *m, int *n);

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

// W = alpha*A*X + beta*Y + gamma
void CSR_MultiplyWithDenseMatrix(
  double *W, int k, int wRowStride, int wColumnStride,
  double alpha, const CSR_Handle *A,
  const double *X, int xRowStride, int xColumnStride,
  double beta, const double *Y, int yRowStride, int yColumnStride,
  double gamma);

void CSR_GetRCMPermutation(const CSR_Handle *A, int *perm, int *inversePerm);
void CSR_GetRCMPermutationWithoutPseudoDiameterSourceSelection(const CSR_Handle *A, int *perm, int *inversePerm);
void CSR_GetBFSPermutation(const CSR_Handle *A, int *perm, int *inversePerm);

void CSR_FindConnectedComponents(
  const CSR_Handle *A,
  int *numOfComponents, int **compToRoot, int **compSizes, int **compSizePrefixSum,
  int **nodesSortedByComp);

void CSR_Permute(const CSR_Handle *A, CSR_Handle *out, const int *columnPerm, const int *rowInversePerm);

int CSR_GetBandwidth(CSR_Handle *A);

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
// sum(x)
double sum(int n, const double *x);
// minimum(x)
double minimum(int n, const double *x);
// w = min(x, alpha)
void min(int n, double *w, const double *x, double alpha); 
void CSR_abs(int n, double *w, const double *x);
void CSR_exp(int n, double *w, const double *x);
void CSR_log1p(int n, double *w, const double *x);

void tallSkinnyDGEMM(
  int transA, int transB,
  int m, int n, int k,   
  double alpha, const double *A, int lda,
  const double *B, int ldb,
  double beta, double *C, int ldc);

#ifdef __cplusplus
}
#endif

#endif // CSR_INTERFACE_H
