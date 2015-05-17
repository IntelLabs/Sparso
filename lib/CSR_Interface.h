#ifndef CSR_INTERFACE_H
#define CSR_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CSR_Handle CSR_Handle;

// const/destruct
CSR_Handle *CSR_Create(int numRows, int numCols, int *i, int *j, double *v, int base);
void CSR_Destroy(CSR_Handle *A);

// accessors
int CSR_GetNumRows(CSR_Handle *A);
int CSR_GetNumCols(CSR_Handle *A);
int CSR_GetNumNonZeros(CSR_Handle *A);

// w = alpha*A*x + beta*y + gamma
void CSR_MultiplyWithVector(
  double *w,
  double alpha, const CSR_Handle *A, const double *x,
  double beta, const double *y,
  double gamma);

void CSR_GetRCMPermutation(const CSR_Handle *A, int *perm, int *inversePerm);
#ifdef USE_BOOST
void CSR_BoostGetRCMPermutation(const CSR_Handle *A, int *perm, int *inversePerm);
void CSR_BoostGetRCMPermutationWithSource(const CSR_Handle *A, int *perm, int *inversePerm, int source);
#endif

void CSR_FindConnectedComponents(
  const CSR_Handle *A,
  int *numOfComponents, int **compToRoot, int **compSizes, int **compSizePrefixSum,
  int **nodesSortedByComp);

void CSR_Permute(const CSR_Handle *A, CSR_Handle *out, const int *columnPerm, const int *rowInversePerm);

int CSR_GetBandwidth(CSR_Handle *A);
void CSR_PrintInDense(CSR_Handle *A);
void CSR_PrintSomeValues(int numRows, int numCols, int *i, int *j, double *v, int distance, bool is_1_based);

void CSR_ReorderMatrix(int numRows, int numCols, int *i, int *j, double *v, int *i1, int *j1, double *v1, 
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
