#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <algorithm>
#include <omp.h>

#include "CSR_Interface.h"
#include "SpMP/mm_io.h"
#include "SpMP/Utils.hpp"

using namespace SpMP;

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
int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: CSR_test matrix_market_file\n");
    return -1;
  }

  {
    int m, n, nnz;
    double *a;
    int *aj, *ai;
    int is_symmetric;
    bool one_based_CSR = false;
    load_matrix_market(argv[1], &ai, &aj, &a, &is_symmetric, &m, &n, &nnz, one_based_CSR, true /*force-symmetric*/);
    printf("m = %d, n = %d, nnz = %d, %csymmetric\n", m, n, nnz, is_symmetric ? ' ' : 'a');
    double bytes = (double)nnz*4;

    CSR_Handle *A = CSR_Create(m, n, ai, aj, a);
    printf("CSR matrix content:\n");
    unsigned int distance = nnz / 100; // print out about 100 elements for manual verification
    //CSR_PrintSomeValues(m, n, ai, aj, a, distance, one_based_CSR);
    printf("original bandwidth: %d\n", CSR_GetBandwidth(A));

    double *x = (double *)malloc(sizeof(double)*n);
    double *y = (double *)malloc(sizeof(double)*m);

    const int REPEAT = 128;

    double t = -omp_get_wtime();
    for (int i = 0; i < REPEAT; ++i) {
      CSR_MultiplyWithVector(y, 1, A, x, 0, y, 0);
    }
    t += omp_get_wtime();

    printf("SpMV BW = %g GB/s\n", ((double)nnz*12 + (m + n)*8)/(t/REPEAT)/1e9);

    int *perm = (int *)malloc(sizeof(int)*m);
    int *inversePerm = (int *)malloc(sizeof(int)*m);
    int source = 1;

    double *a2 = (double *)malloc(sizeof(double)*nnz);
    int *aj2 = (int *)malloc(sizeof(int)*nnz);
    int *ai2 = (int *)malloc(sizeof(int)*(m + 1));
    CSR_Handle *A2 = CSR_Create(m, n, ai2, aj2, a2);

    for (int permuteType = 0; permuteType < 3; ++permuteType) {
      if (2 == permuteType) {
        printf("RCM permutation\n");
      }
      else if (1 == permuteType) {
        printf("RCM permutation w/o source selection\n");
      }
      else {
        printf("BFS permutation\n");
      }

      t = -omp_get_wtime();
      if (2 == permuteType) {
        CSR_GetRCMPermutation(A, perm, inversePerm);
      }
      else if (1 == permuteType) {
        CSR_GetRCMPermutationWithoutPseudoDiameterSourceSelection(A, perm, inversePerm);
      }
      else {
        CSR_GetBFSPermutation(A, perm, inversePerm);
      }
      t += omp_get_wtime();
      printf("Constructing permutation takes %f (%f GB/s)\n", t, bytes/t/1e9);

      SpMP::isPerm(perm, m);
      SpMP::isPerm(inversePerm, m);

      t = -omp_get_wtime();
      CSR_Permute(A, A2, perm, inversePerm);
      t += omp_get_wtime();
      printf("Permute takes %f (%f GB/s)\n", t, bytes/t/1e9);
      printf("Permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

      t = -omp_get_wtime();
      for (int i = 0; i < REPEAT; ++i) {
        CSR_MultiplyWithVector(y, 1, A2, x, 0, y, 0);
      }
      t += omp_get_wtime();
      printf("SpMV BW = %g GB/s\n", ((double)nnz*12 + (m + n)*8)/(t/REPEAT)/1e9);
    }

    CSR_Destroy(A2);
    CSR_Destroy(A);
  }

#if 0
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
    CSR_GetRCMPermutation(A, perm, inversePerm);
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
    CSR_GetRCMPermutation(A, perm, inversePerm);
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
#endif

  return 0;
}
