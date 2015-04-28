#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <algorithm>
#include <omp.h>

#include "CSR_Interface.h"
#include "mm_io.h"

bool isPerm(int *perm, int n)
{
  int *temp = new int[n];
  memcpy(temp, perm, sizeof(int)*n);
  std::sort(temp, temp + n);
  int *last = std::unique(temp, temp + n);
  if (last != temp + n) {
    memcpy(temp, perm, sizeof(int)*n);
    std::sort(temp, temp + n);

    for (int i = 0; i < n; ++i) {
      if (temp[i] == i - 1) {
        printf("%d duplicated\n", i);
        assert(false);
        return false;
      }
      else if (temp[i] != i) {
        printf("%d missed\n", i);
        assert(false);
        return false;
      }
    }
  }
  delete[] temp;
  return true;
}

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
    load_matrix_market(argv[1], &a, &aj, &ai, &is_symmetric, &m, &n, &nnz, one_based_CSR);
    printf("m = %d, n = %d, nnz = %d, %csymmetric\n", m, n, nnz, is_symmetric ? ' ' : 'a');
    double bytes = (double)nnz*12;

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
      CSR_MultiplyWithVector(A, y, x);
    }
    t += omp_get_wtime();

    printf("SpMV BW = %g GB/s\n", ((double)nnz*12 + (m + n)*8)/(t/REPEAT)/1e9);

    int *perm = (int *)malloc(sizeof(int)*m);
    int *inversePerm = (int *)malloc(sizeof(int)*m);
    int source = 1;

    double *a2 = (double *)malloc(sizeof(double)*nnz);
    int *aj2 = (int *)malloc(sizeof(int)*nnz);
    int *ai2 = (int *)malloc(sizeof(int)*(m + 1));
    CSR_Handle *A2 = NULL;

    printf("RCM permutation\n");

#ifdef USE_BOOST
    t = -omp_get_wtime();
    CSR_BoostGetRCMPemutation(A, perm, inversePerm);
    t += omp_get_wtime();
    isPerm(perm, m);
    isPerm(inversePerm, m);

    printf("Boost RCM takes %f (%f GB/s)\n", t, bytes/t/1e9);

    A2 = CSR_Create(m, n, ai2, aj2, a2);

    CSR_Permute(A, A2, perm, inversePerm);
    printf("Boost RCM permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));
#endif

    /*t = -omp_get_wtime();
    CSR_GetRCMPemutationWithSource(A, perm, inversePerm, source);
    t += omp_get_wtime();
    isPerm(perm, m);
    isPerm(inversePerm, m);
    
    printf("Boost RCM with source %d takes %f (%f GB/s)\n", t, source, bytes/t/1e9);

    CSR_Permute(A, A2, perm, inversePerm);
    printf("RCM permuted bandwidth with source %d: %d\n\n", source, CSR_GetBandwidth(A2));

    t = -omp_get_wtime();
    CSR_GetRCMPemutationNewWithSource(A, perm, inversePerm, source);
    t += omp_get_wtime();
    isPerm(perm, m);
    isPerm(inversePerm, m);
    
    printf("My RCM with source %d takes %f (%f GB/s)\n", source, t, bytes/t/1e9);

    CSR_Permute(A, A2, perm, inversePerm);
    printf("My RCM permuted bandwidth with source %d: %d\n\n", source, CSR_GetBandwidth(A2));*/

    t = -omp_get_wtime();
    CSR_GetRCMPemutation(A, perm, inversePerm);
    t += omp_get_wtime();
    isPerm(perm, m);
    isPerm(inversePerm, m);
    
    printf("My RCM takes %f (%f GB/s)\n", t, bytes/t/1e9);

    A2 = CSR_Create(m, n, ai2, aj2, a2);

    CSR_Permute(A, A2, perm, inversePerm);
    printf("My permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

    t = -omp_get_wtime();
    for (int i = 0; i < REPEAT; ++i) {
      CSR_MultiplyWithVector(A2, y, x);
    }
    t += omp_get_wtime();

    printf("SpMV BW = %g GB/s\n", ((double)nnz*12 + (m + n)*8)/(t/REPEAT)/1e9);

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
#endif

  return 0;
}
