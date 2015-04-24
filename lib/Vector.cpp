#include <cstdio>
#include <cassert>
#include <algorithm>

#include <omp.h>

using namespace std;

extern "C" void reorderVector(double *v, double *tmp, const int *perm, int len)
{
  if (!perm) return;

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    tmp[perm[i]] = v[i];
  }

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    v[i] = tmp[i];
  }
}

extern "C" void reorderVectorWithInversePerm(double *v, double *tmp, const int *inversePerm, int len)
{
  if (!inversePerm) return;

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    tmp[i] = v[inversePerm[i]];
  }

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    v[i] = tmp[i];
  }
}

  /**
   * Compute y = A*x
   */
// TODO: remove this once MKL libray call is fine, or when reusing 
// works so that we can convert 0 to 1 based only once in the loop
// This is a temporary workaround. To remove in future.
extern "C" void CSR_MultiplyWithVector_1Based(int num_rows, int *rowPtr, int *colIdx, double* values, double *x, double *y)
{
//#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
#endif

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    int nnz = rowPtr[num_rows] - 1;
    int nnzPerThread = (nnz + nthreads - 1)/nthreads;
    int iBegin = lower_bound(rowPtr, rowPtr + num_rows, nnzPerThread*tid + 1) - rowPtr;
    int iEnd = lower_bound(rowPtr, rowPtr + num_rows, nnzPerThread*(tid + 1) + 1) - rowPtr;
    assert(iBegin <= iEnd);
    assert(iBegin >= 0 && iBegin <= num_rows);
    assert(iEnd >= 0 && iEnd <= num_rows);

    for (int i = iBegin; i < iEnd; ++i) {
      double sum = 0;
      for (int j = rowPtr[i] - 1; j < rowPtr[i + 1] - 1; ++j) {
        sum += values[j]*x[colIdx[j] - 1];
      }
      y[i] = sum;
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
#endif
  } // omp parallel
}
