#include <stdio.h>
#include <assert.h>

#include "CSR_Interface.h"

/* Expected output

original banwdith: 8
RCM permutation
0 8 5 7 3 6 4 2 1 9
RCM permuted bandwidth: 4

original banwdith: 6
RCM permutation
3 6 0 4 2 1 7 5
RCM permuted bandwidth: 2
 */
int main()
{
  {
    // from http://www.boost.org/doc/libs/1_37_0/libs/graph/example/cuthill_mckee_ordering.cpp
    int colIdx[] =    { 3, 5, 2, 4, 6, 9, 1, 3, 4, 0, 2, 5, 8, 1, 2, 6, 0, 3, 6, 7, 1, 4, 5, 7, 5, 6, 3, 1 };
    int rowPtr[] =    { 0,    2,          6,       9,          13,      16,         20,         24,   26,27, 28 };
    double values[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    int m = sizeof(rowPtr)/sizeof(rowPtr[0]) - 1;
    assert(10 == m);
    int nnz = rowPtr[m];
    assert(nnz == sizeof(colIdx)/sizeof(colIdx[0]));
    CSR_Handle *A = CSR_Create(m, m, rowPtr, colIdx, values);
    CSR_PrintInDense(A);
    printf("original banwdith: %d\n", CSR_GetBandwidth(A));

    printf("RCM permutation\n");
    int perm[m], inversePerm[m];
    CSR_GetRCMPemutation(A, perm, inversePerm);
    for (int i = 0; i < m; ++i) {
      printf("%d ", inversePerm[i]);
    }
    printf("\n");

    int rowPtr2[m + 1];
    int colIdx2[nnz];
    double values2[nnz];
    CSR_Handle *A2 = CSR_Create(m, m, rowPtr2, colIdx2, values2);

    CSR_Permute(A, A2, perm, inversePerm);
    CSR_PrintInDense(A2);
    printf("RCM permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

    CSR_Destroy(A2);
    CSR_Destroy(A);
  }

  {
    // from http://ciprian-zavoianu.blogspot.com/2009/01/project-bandwidth-reduction.html
    int colIdx[] =    { 0, 4, 1, 2, 5, 7, 1, 2, 4, 3, 6, 0, 2, 4, 1, 5, 7, 3, 6, 1, 5, 7 };
    int rowPtr[] =    { 0,    2,          6,       9,    11,      14,      17,   19,     22 };
    double values[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    int m = sizeof(rowPtr)/sizeof(rowPtr[0]) - 1;
    int nnz = rowPtr[m];
    assert(nnz == sizeof(colIdx)/sizeof(colIdx[0]));
    CSR_Handle *A = CSR_Create(m, m, rowPtr, colIdx, values);
    CSR_PrintInDense(A);
    printf("original banwdith: %d\n", CSR_GetBandwidth(A));

    printf("RCM permutation\n");
    int perm[m], inversePerm[m];
    CSR_GetRCMPemutation(A, perm, inversePerm);
    for (int i = 0; i < m; ++i) {
      printf("%d ", inversePerm[i]);
    }
    printf("\n");

    int rowPtr2[m + 1];
    int colIdx2[nnz];
    double values2[nnz];
    CSR_Handle *A2 = CSR_Create(m, m, rowPtr2, colIdx2, values2);

    CSR_Permute(A, A2, perm, inversePerm);
    CSR_PrintInDense(A2);
    printf("RCM permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

    CSR_Destroy(A2);
    CSR_Destroy(A);
  }

  return 0;
}

