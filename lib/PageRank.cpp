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


