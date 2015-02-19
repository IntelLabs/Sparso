#include "CSR.hpp"
#include "CSR_Interface.h"

extern "C" {

CSR_Handle *CSR_Create(int numRows, int numCols, int *i, int *j, double *v)
{
  return (CSR_Handle *)(new CSR(numRows, numCols, i, j, v));
}

void CSR_MultiplyWithVector(const CSR_Handle *A, double *y, const double *x)
{
  ((CSR *)A)->multiplyWithVector(y, x);
}

void CSR_GetRCMPemutation(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm);
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

}
