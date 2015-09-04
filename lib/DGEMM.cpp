#include <mkl.h>

#include "SpMP/Utils.hpp"

#include "CSR_Interface.h"

using namespace SpMP;

extern "C" {

// all matrices are assumed to be in column major
void tallSkinnyDGEMM(
  int transA, int transB,
  int m, int n, int k,   
  double alpha, const double *A, int lda,
  const double *B, int ldb,
  double beta, double *C, int ldc)
{ 
  if (!transA && !transB && m > n && m > k) {
    // tall skinny matrix multiplied by small matrix
#pragma omp parallel
    {
      int begin, end;
      getSimpleThreadPartition(&begin, &end, m);

      cblas_dgemm(
        CblasColMajor, CblasNoTrans, transB ? CblasTrans : CblasNoTrans,
        end - begin, n, k,
        alpha, A + begin, lda,
        B, ldb,
        beta, C + begin, ldc);
    }   
  } 
  else if (!transB && m <= 24 && n <= 24) {
    double *temp = MALLOC(double, m*n*(omp_get_max_threads() - 1));

#pragma omp parallel
    {
      int begin, end;
      getSimpleThreadPartition(&begin, &end, k);

      int tid = omp_get_thread_num();

      cblas_dgemm(
        CblasColMajor,
        transA ? CblasTrans : CblasNoTrans,
        transB ? CblasTrans : CblasNoTrans,
        m, n, end - begin,
        alpha, A + (transA ? begin : begin*lda), lda,
        B + begin, ldb,  
        0 == tid ? beta : 0,
        0 == tid ? C : temp + m*n*(tid - 1),
        0 == tid ? ldc : m);
    }

    for (int t = 0; t < omp_get_max_threads() - 1; ++t) {
      for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
          C[i + j*ldc] += temp[i + j*m + t*m*n];
        }
      }
    } // for each thread

    FREE(temp)
  }
  else {
    cblas_dgemm(
      CblasColMajor,
      transA ? CblasTrans : CblasNoTrans,
      transB ? CblasTrans : CblasNoTrans,
      m, n, k,
      alpha, A, lda,
      B, ldb,
      beta, C, ldc);
  }
}

} // extern "C"
