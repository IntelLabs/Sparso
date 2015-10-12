#include <algorithm>

#include "SpMP/CSR.hpp"

#include "CSR_Interface.h"

using namespace std;
using namespace SpMP;

// w = (alpha*A*x + beta*y + gamma)*z
template<class T, bool HAS_VALUE = true, bool HAS_Z = false>
inline
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
          if (HAS_Z) {
            if (HAS_VALUE) {
              for (int i = iBegin; i < iEnd; ++i) {
                T sum = 0;
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
                  sum += values[j]*x[colidx[j]];
                }
                w[i] = (alpha*sum + gamma)*z[i];
              }
            }
            else {
              for (int i = iBegin; i < iEnd; ++i) {
                T sum = 0;
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
                  sum += x[colidx[j]];
                }
                w[i] = (alpha*sum + gamma)*z[i];
              }
            }
          }
          else {
            if (HAS_VALUE) {
              for (int i = iBegin; i < iEnd; ++i) {
                T sum = 0;
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
                  sum += values[j]*x[colidx[j]];
                }
                w[i] = alpha*sum + gamma;
              }
            }
            else {
              for (int i = iBegin; i < iEnd; ++i) {
                T sum = 0;
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
                  sum += x[colidx[j]];
                }
                w[i] = alpha*sum + gamma;
              }
            }
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

extern double *getTempVector(int l); // defined in knob.cpp

extern "C" {

// context insensitive version of SpMV
// Since we don't know if A in CSC is symmetric
// we resort into SpMV in CSC.
// transpose to CSR and then performing CSR SpMV
// turns out to be slower
void CSR_MultiplyWithVector(
  double *w,
  double alpha, const CSR_Handle *A, const double *x,
  double beta, const double *y,
  double gamma,
  const double *z)
{
  CSR *AT = ((CSR *)A);

  int base = AT->getBase();

  int *rowptr = AT->rowptr - base;
  int *colidx = AT->colidx - base;
  double *values = AT->values - base;

  x -= base;

//#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
#endif
  
  double *temp_buffer_array = MALLOC(double, omp_get_max_threads()*AT->n);

#pragma omp parallel
  {
    int iBegin, iEnd;
    getLoadBalancedPartition(&iBegin, &iEnd, rowptr + base, AT->m);
    iBegin += base;
    iEnd += base;

    int tid = omp_get_thread_num();
    double *temp_buffer = temp_buffer_array + tid*AT->n;
    for (int i = 0; i < AT->n; ++i) {
      temp_buffer[i] = 0;
    }
    temp_buffer -= base;

    for (int i = iBegin; i < iEnd; ++i) {
      double xi = x[i];
      for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
        temp_buffer[colidx[j]] += xi*values[j];
      }
    }

#pragma omp barrier

#pragma omp for
    for (int i = 0; i < AT->n; ++i) {
      double sum = temp_buffer_array[i];
      for (int j = 1; j < omp_get_num_threads(); ++j) {
        sum += temp_buffer_array[j*AT->n + i];
      }
      if (z) {
        if (beta == 0) {
          w[i] = (alpha*sum + gamma)*z[i];
        }
        else {
          w[i] = (alpha*sum + beta*y[i] + gamma)*z[i];
        }
      }
      else {
        if (beta == 0) {
          w[i] = alpha*sum + gamma;
        }
        else {
          w[i] = alpha*sum + beta*y[i] + gamma;
        }
      }
    }
  }

  FREE(temp_buffer_array);
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
