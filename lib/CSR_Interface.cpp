#include "CSR.hpp"
#include "CSR_Interface.h"

extern "C" {

CSR_Handle *CSR_Create(int numRows, int numCols, int *i, int *j, double *v)
{
  return (CSR_Handle *)(new CSR(numRows, numCols, i, j, v));
}

int CSR_GetNumRows(CSR_Handle *A)
{
  return ((CSR *)A)->m;
}

int CSR_GetNumCols(CSR_Handle *A)
{
  return ((CSR *)A)->n;
}

int CSR_GetNumNonZeros(CSR_Handle *A)
{
  return ((CSR *)A)->rowPtr[((CSR *)A)->m];
}

void CSR_MultiplyWithVector(const CSR_Handle *A, double *y, const double *x)
{
  ((CSR *)A)->multiplyWithVector(y, x);
}

void CSR_GetRCMPemutation(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm);
}

void CSR_GetRCMPemutationWithSource(const CSR_Handle *A, int *perm, int *inversePerm, int source)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm, source);
}

void CSR_GetRCMPemutationNew(const CSR_Handle *A, int *perm, int *inversePerm, int source)
{
  ((CSR *)A)->getRCMPermutationNew(perm, inversePerm, source);
}

void CSR_Permute(const CSR_Handle *A, CSR_Handle *out, const int *columnPerm, const int *rowInversePerm)
{
  ((CSR *)A)->permute((CSR *)out, columnPerm, rowInversePerm);
}

int CSR_GetBandwidth(CSR_Handle *A)
{
  return ((CSR *)A)->getBandwidth();
}

void CSR_PrintInDense(CSR_Handle *A)
{
  ((CSR *)A)->printInDense();
}

void CSR_Destroy(CSR_Handle *A)
{
  delete (CSR *)A;
}

// The first few parameters (numRows to v) represent the source matrix to be reordered. The next few
// parameters (i1~v1) are the spaces that have been allocated to store the results. 
// Perm and inversePerm are the spaces that have been allocated for permutation and inverse permutation
// info; when getPermutation is true, this function computes and stores the info into them; otherwise,
// they already contain the info, and this function just uses it.
void CSR_ReorderMatrix(int numRows, int numCols, int *i, int *j, double *v, int *i1, int *j1, double *v1, 
                 int *perm, int *inversePerm, bool getPermutation)
{
    CSR_Handle *A = CSR_Create(numRows, numCols, i, j, v);
    if (getPermutation) {
        CSR_GetRCMPemutation(A, perm, inversePerm);
    }
    CSR_Handle *newA = CSR_Create(numRows, numCols, i1, j1, v1);
    CSR_Permute(A, newA, perm, inversePerm);
    CSR_Destroy(newA);
    CSR_Destroy(A);
}

}
