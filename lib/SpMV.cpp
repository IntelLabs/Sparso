#include <algorithm>

#include "SpMP/CSR.hpp"

#include "CSR_Interface.h"

using namespace std;
using namespace SpMP;

// w = (alpha*A*x + beta*y + gamma)*z
template<class T, bool HAS_VALUE = true, bool HAS_Z = false>
static void SpMV_(
  int m,
  T *w,
  T alpha,
  const int *rowptr, const int *colidx, const T* values,
  const T *x,
  T beta,
  const T *y,
  T gamma,
  const T *z)
{
  assert(w != x);

  int base = rowptr[0];

  rowptr -= base;
  colidx -= base;
  if (values) values -= base;

  w -= base;
  x -= base;
  if (y) y -= base;
  if (z) z -= base;

//#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
#endif

#pragma omp parallel
  {
    int iBegin, iEnd;
    getLoadBalancedPartition(&iBegin, &iEnd, rowptr + base, m);
    iBegin += base;
    iEnd += base;

    if (1 == alpha) {
      if (0 == beta) {
        if (0 == gamma) {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
        else {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum += gamma;
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
      }
      else {
        // alpha == 1 && beta != 0
        if (0 == gamma) {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum += beta*y[i];
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
        else {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum += beta*y[i] + gamma;
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
      }
    }
    else {
      // alpha != 1
      if (0 == beta) {
        if (0 == gamma) {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum *= alpha;
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
        else {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum = alpha*sum + beta*y[i];
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
      }
      else {
        // alpha != 1 && beta != 0
        if (0 == gamma) {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum = alpha*sum + beta*y[i];
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
        else {
          for (int i = iBegin; i < iEnd; ++i) {
            T sum = 0;
            for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
              if (HAS_VALUE) {
                sum += values[j]*x[colidx[j]];
              }
              else {
                sum += x[colidx[j]];
              }
            }
            sum = alpha*sum + beta*y[i] + gamma;
            if (HAS_Z) w[i] = z[i]*sum;
            else w[i] = sum;
          }
        }
      }
    } // alpha != 1

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
}

void multiplyWithVector(
  double *w,
  double alpha, const CSR *A, const double *x,
  double beta, const double *y,
  double gamma,
  const double *z)
{
  if (A->values) {
    if (z) {
      SpMV_<double, true, true>(
        A->m, w, alpha, A->rowptr, A->colidx, A->values, x, beta, y, gamma, z);
    }
    else {
      SpMV_<double, true, false>(
        A->m, w, alpha, A->rowptr, A->colidx, A->values, x, beta, y, gamma, z);
    }
  }
  else {
    if (z) {
      SpMV_<double, false, true>(
        A->m, w, alpha, A->rowptr, A->colidx, A->values, x, beta, y, gamma, z);
    }
    else {
      SpMV_<double, false, false>(
        A->m, w, alpha, A->rowptr, A->colidx, A->values, x, beta, y, gamma, z);
    }
  }
}

extern "C" {

void CSR_MultiplyWithVector(
  double *w,
  double alpha, const CSR_Handle *A, const double *x,
  double beta, const double *y,
  double gamma,
  const double *z)
{
  multiplyWithVector(w, alpha, (CSR *)A, x, beta, y, gamma, z);
}

void CSR_MultiplyWithDenseMatrix(
  double *W, int k, int wRowStride, int wColumnStride,
  double alpha, const CSR_Handle *A,
  const double *X, int xRowStride, int xColumnStride,
  double beta, const double *Y, int yRowStride, int yColumnStride,
  double gamma)
{
  ((CSR *)A)->multiplyWithDenseMatrix(
    W, k, wRowStride, wColumnStride,
    alpha,
    X, xRowStride, xColumnStride,
    beta, Y, yRowStride, yColumnStride,
    gamma);
}

} // extern "C"
