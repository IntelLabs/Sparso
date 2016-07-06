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
#include <cassert>
#include <algorithm>

#include <omp.h>

#include "CSR_Interface.h"

using namespace std;

extern "C" {
// w = (alpha*A*x + beta*y + gamma).*z
void PageRank(
  int num_rows,
  double *w,
  double alpha,
  const int *rowPtr, const int *colIdx,
  const double *x,
  double beta,
  const double *y,
  double gamma,
  const double *z)
{
  const int BASE = 1;

//#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
#endif

#ifdef SEP
  static int cnt = 0;
  if (2 == cnt) {
    VTResumeSampling();
  }
#endif

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    int nnz = rowPtr[num_rows] - BASE;
    int nnzPerThread = (nnz + nthreads - 1)/nthreads;
    int iBegin = lower_bound(rowPtr, rowPtr + num_rows, nnzPerThread*tid + BASE) - rowPtr;
    int iEnd = lower_bound(rowPtr, rowPtr + num_rows, nnzPerThread*(tid + 1) + BASE) - rowPtr;
    assert(iBegin <= iEnd);
    assert(iBegin >= 0 && iBegin <= num_rows);
    assert(iEnd >= 0 && iEnd <= num_rows);

    for (int i = iBegin; i < iEnd; ++i) {
      double sum = 0;
      for (int j = rowPtr[i] - BASE; j < rowPtr[i + 1] - BASE; ++j) {
        sum += x[colIdx[j] - BASE];
      }
      w[i] = (alpha*sum + beta*y[i] + gamma)*z[i];
    }

#ifdef MEASURE_LOAD_BALANCE
    double t = omp_get_wtime();
#pragma omp barrier
    barrierTimes[tid] = omp_get_wtime() - t;

#pragma omp barrier
#pragma omp master
    {
      double tEnd = omp_get_wtime();
      double barrierTimeSum = 0;
      for (int i = 0; i < nthreads; ++i) {
        barrierTimeSum += barrierTimes[i];
      }
      printf("%f load imbalance = %f\n", tEnd - tBegin, barrierTimeSum/(tEnd - tBegin)/nthreads);
    }
#undef MEASURE_LOAD_BALANCE
#endif // MEASURE_LOAD_BALANCE
  } // omp parallel

#ifdef SEP
  if (2 == cnt) {
    VTPauseSampling();
  }
  ++cnt;
#endif
}
} // extern "C"


