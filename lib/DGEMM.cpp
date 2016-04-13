/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <mkl.h>
#include <immintrin.h>

#include "SpMP/Utils.hpp"

#include "CSR_Interface.h"

using namespace SpMP;

template<int n, int k, int alpha, int beta>
static void tallSkinnyDGEMM1_(
  int m,
  const double *A, int lda,
  const double *B, int ldb,
  double *C, int ldc)
{

  // tall skinny matrix multiplied by small matrix
#pragma omp parallel
  {
    const int IB = 32;

    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    int nblocks = (m + IB - 1)/IB;
    int blocksPerThread = (nblocks + nthreads - 1)/nthreads;

    int iBegin = std::min(blocksPerThread*IB*tid, m);
    int iEnd = std::min(iBegin + blocksPerThread*IB, m);
    int iEndFloor = iBegin + (iEnd - iBegin)/IB*IB;

    for (int i = iBegin; i < iEndFloor; i += IB) {
//#pragma unroll_and_jam (n)
      for (int j = 0; j < n; j++) {
        __m256d sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7;
        if (beta == 0) {
          sum0 = _mm256_setzero_pd();
          sum1 = _mm256_setzero_pd();
          sum2 = _mm256_setzero_pd();
          sum3 = _mm256_setzero_pd();
          sum4 = _mm256_setzero_pd();
          sum5 = _mm256_setzero_pd();
          sum6 = _mm256_setzero_pd();
          sum7 = _mm256_setzero_pd();
        }
        else if (alpha == -1) {
          if (beta == 1) {
            sum0 = _mm256_load_pd(C + i + j*ldc);
            sum1 = _mm256_load_pd(C + i + 4 + j*ldc);
            sum2 = _mm256_load_pd(C + i + 8 + j*ldc);
            sum3 = _mm256_load_pd(C + i + 12 + j*ldc);
            sum4 = _mm256_load_pd(C + i + 16 + j*ldc);
            sum5 = _mm256_load_pd(C + i + 20 + j*ldc);
            sum6 = _mm256_load_pd(C + i + 24 + j*ldc);
            sum7 = _mm256_load_pd(C + i + 28 + j*ldc);
          }
          else {
            sum0 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + j*ldc));
            sum1 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 4 + j*ldc));
            sum2 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 8 + j*ldc));
            sum3 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 12 + j*ldc));
            sum4 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 16 + j*ldc));
            sum5 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 20 + j*ldc));
            sum6 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 24 + j*ldc));
            sum7 = _mm256_mul_pd(
              _mm256_set1_pd(beta),
              _mm256_load_pd(C + i + 28 + j*ldc));
          }
        }
        else {
          sum0 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + j*ldc));
          sum1 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 4 + j*ldc));
          sum2 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 8 + j*ldc));
          sum3 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 12 + j*ldc));
          sum4 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 16 + j*ldc));
          sum5 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 20 + j*ldc));
          sum6 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 24 + j*ldc));
          sum7 = _mm256_mul_pd(
            _mm256_set1_pd(beta/alpha),
            _mm256_load_pd(C + i + 28 + j*ldc));
        }

#pragma unroll (k)
        for (int l = 0; l < k; l++) {
          __m256d b_lj = _mm256_set1_pd(B[l + j*ldb]);

          if (-1 == alpha) {
            sum0 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + l*lda), b_lj, sum0);
            sum1 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 4 + l*lda), b_lj, sum1);
            sum2 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 8 + l*lda), b_lj, sum2);
            sum3 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 12 + l*lda), b_lj, sum3);
            sum4 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 16 + l*lda), b_lj, sum4);
            sum5 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 20 + l*lda), b_lj, sum5);
            sum6 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 24 + l*lda), b_lj, sum6);
            sum7 = _mm256_fnmadd_pd(_mm256_load_pd(A + i + 28 + l*lda), b_lj, sum7);
          }
          else {
            sum0 = _mm256_fmadd_pd(_mm256_load_pd(A + i + l*lda), b_lj, sum0);
            sum1 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 4 + l*lda), b_lj, sum1);
            sum2 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 8 + l*lda), b_lj, sum2);
            sum3 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 12 + l*lda), b_lj, sum3);
            sum4 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 16 + l*lda), b_lj, sum4);
            sum5 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 20 + l*lda), b_lj, sum5);
            sum6 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 24 + l*lda), b_lj, sum6);
            sum7 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 28 + l*lda), b_lj, sum7);
          }
        }

        if (1 == alpha || -1 == alpha) {
          _mm256_store_pd(C + i + j*ldc, sum0);
          _mm256_store_pd(C + i + 4 + j*ldc, sum1);
          _mm256_store_pd(C + i + 8 + j*ldc, sum2);
          _mm256_store_pd(C + i + 12 + j*ldc, sum3);
          _mm256_store_pd(C + i + 16 + j*ldc, sum4);
          _mm256_store_pd(C + i + 20 + j*ldc, sum5);
          _mm256_store_pd(C + i + 24 + j*ldc, sum6);
          _mm256_store_pd(C + i + 28 + j*ldc, sum7);
        }
        else {
          _mm256_store_pd(
            C + i + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum0));
          _mm256_store_pd(
            C + i + 4 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum1));
          _mm256_store_pd(
            C + i + 8 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum2));
          _mm256_store_pd(
            C + i + 12 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum3));
          _mm256_store_pd(
            C + i + 16 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum4));
          _mm256_store_pd(
            C + i + 20 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum5));
          _mm256_store_pd(
            C + i + 24 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum6));
          _mm256_store_pd(
            C + i + 28 + j*ldc,
            _mm256_mul_pd(_mm256_set1_pd(alpha), sum7));
        }
      }
    }

    for (int i = iEndFloor; i < iEnd; ++i) {
      for (int j = 0; j < n; j++) {
        double sum = 0;
        for (int l = 0; l < k; l++) {
          sum += A[i + l*lda]*B[l + j*ldb];
        }
        C[i + j*ldc] = beta*C[i + j*ldc] + alpha*sum;
      }
    }
  } // omp parallel
}

extern "C" {

unsigned long long t1 = 0, t2 = 0, t3 = 0;

void printStat() {
  printf("%f %f %f\n", t1/get_cpu_freq(), t2/get_cpu_freq(), t3/get_cpu_freq());
}

// all matrices are assumed to be in column major
void tallSkinnyDGEMM(
  int transA, int transB,
  int m, int n, int k,   
  double alpha, const double *A, int lda,
  const double *B, int ldb,
  double beta, double *C, int ldc)
{ 
  unsigned long long t_begin = __rdtsc();

  if (!transA && !transB && m > n && m > k) {
    // tall skinny matrix multiplied by small matrix

    bool executed = false;
    if (n == 3) {
      if (k == 1 && alpha == -1 && beta == 1) {
        tallSkinnyDGEMM1_<3, 1, -1, 1>(m, A, lda, B, ldb, C, ldc);
        executed = true;
      }
      else if (k == 2 && alpha == -1 && beta == 1) {
        tallSkinnyDGEMM1_<3, 2, -1, 1>(m, A, lda, B, ldb, C, ldc);
        executed = true;
      }
      else if (k == 3 && alpha == -1 && beta == 1) {
        tallSkinnyDGEMM1_<3, 3, -1, 1>(m, A, lda, B, ldb, C, ldc);
        executed = true;
      }
      else if (k == 6 && alpha == -1 && beta == 1) {
        tallSkinnyDGEMM1_<3, 6, -1, 1>(m, A, lda, B, ldb, C, ldc);
        executed = true;
      }
      else if (k == 9 && alpha == -1 && beta == 1) {
        tallSkinnyDGEMM1_<3, 9, -1, 1>(m, A, lda, B, ldb, C, ldc);
        executed = true;
      }
    }
    else if (n == 9) {
      if (k == 3 && alpha == 1 && beta == 0) {
        tallSkinnyDGEMM1_<9, 3, 1, 0>(m, A, lda, B, ldb, C, ldc);
        executed = true;
      }
    }

    if (!executed) {
      printf("1 %f %f %d %d\n", alpha, beta, n, k);
#pragma omp parallel
      {
        const int IB = 8;

        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        int nblocks = (m + IB - 1)/IB;
        int blocksPerThread = (nblocks + nthreads - 1)/nthreads;

        int iBegin = std::min(blocksPerThread*IB*tid, m);
        int iEnd = std::min(iBegin + blocksPerThread*IB, m);
        int iEndFloor = iBegin + (iEnd - iBegin)/IB*IB;

        cblas_dgemm(
          CblasColMajor, CblasNoTrans, transB ? CblasTrans : CblasNoTrans,
          iEnd - iBegin, n, k,
          alpha, A + iBegin, lda,
          B, ldb,
          beta, C + iBegin, ldc);

        /*for (int i = iBegin; i < iEndFloor; i += IB) {
          for (int j = 0; j < n; j++) {
            __m256d sum0 = _mm256_mul_pd(
              _mm256_set1_pd(beta/alpha),
              _mm256_load_pd(C + i + j*ldc));
            __m256d sum1 = _mm256_mul_pd(
              _mm256_set1_pd(beta/alpha),
              _mm256_load_pd(C + i + 4 + j*ldc));

            for (int l = 0; l < k; l++) {
              sum0 = _mm256_fmadd_pd(_mm256_load_pd(A + i + l*lda), _mm256_set1_pd(B[l + j*ldb]), sum0);
              sum1 = _mm256_fmadd_pd(_mm256_load_pd(A + i + 4 + l*lda), _mm256_set1_pd(B[l + j*ldb]), sum1);
            }

            _mm256_store_pd(
              C + i + j*ldc,
              _mm256_mul_pd(_mm256_set1_pd(alpha), sum0));
            _mm256_store_pd(
              C + i + 4 + j*ldc,
              _mm256_mul_pd(_mm256_set1_pd(alpha), sum1));
          }
        }

        for (int i = iEndFloor; i < iEnd; ++i) {
          for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int l = 0; l < k; l++) {
              sum += A[i + l*lda]*B[l + j*ldb];
            }
            C[i + j*ldc] = beta*C[i + j*ldc] + alpha*sum;
          }
        }*/
      } // omp parallel
    }

    t1 += __rdtsc() - t_begin;
  } 
  else if (!transB && m <= 24 && n <= 24) {
    // fat matrix multiplied by skinny matrix
    if (!transA) {
      printf("2 %d %d %d %d %d\n", transA, transB, m, k, n);
    }
    //double *temp = MALLOC(double, m*n*(omp_get_max_threads() - 1));

#pragma omp parallel for collapse(2)
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        double sum = 0;
        for (int l = 0; l < k; ++l) {
          sum += A[i*lda + l]*B[j*ldb + l];
        }
        C[i + j*ldc] = alpha*sum + beta*C[i + j*ldc];
      }
    }

/*#pragma omp parallel
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

    FREE(temp)*/

    t2 += __rdtsc() - t_begin;
  }
  else {
    //printf("3 %d %d %d\n", m, k, n);
    cblas_dgemm(
      CblasColMajor,
      transA ? CblasTrans : CblasNoTrans,
      transB ? CblasTrans : CblasNoTrans,
      m, n, k,
      alpha, A, lda,
      B, ldb,
      beta, C, ldc);

    t3 += __rdtsc() - t_begin;
  }
}

} // extern "C"
