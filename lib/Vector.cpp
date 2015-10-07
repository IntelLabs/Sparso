#include <cstdio>
#include <cassert>
#include <cfloat>
#include <algorithm>

#include <omp.h>
#include <mkl.h>

#include "SpMP/Utils.hpp"

using namespace std;

// w = alpha*x + beta*w
template<class T>
void waxpby_(int n, T *w, T alpha, const T *x, T beta, const T *y)
{
  if (1 == alpha) {
    if (1 == beta) {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = x[i] + y[i];
      }
    }
    else if (-1 == beta) {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = x[i] - y[i];
      }
    }
    else {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = x[i] + beta*y[i];
      }
    }
  }
  else if (1 == beta) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i] + y[i];
    }
  }
  else if (-1 == alpha) {
    if (0 == beta) {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = -x[i];
      }
    }
    else {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = beta*y[i] - x[i];
      }
    }
  }
  else if (0 == beta) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i];
    }
  }
  else {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i] + beta*y[i];
    }
  }
}

// w = alpha*x + beta
template<class T>
void waxpb_(int n, T *w, T alpha, const T *x, T beta)
{
  if (1 == alpha) {
    if (0 == beta) {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = x[i];
      }
    }
    else {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        w[i] = x[i] + beta;
      }
    }
  }
  else if (0 == beta) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i];
    }
  }
  else {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i] + beta;
    }
  }
}

template<class T>
T dot_(int n, const T *x, const T *y)
{
  T sum = 0;
#pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; ++i) {
    sum += x[i]*y[i];
  }
  return sum;
}

template<class T>
void pointwiseDivide_(int n, T *w, const T *x, const T *y)
{
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    w[i] = x[i]/y[i];
  }
}

template<class T>
void pointwiseMultiply_(int n, T *w, const T *x, const T *y)
{
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    w[i] = x[i]*y[i];
  }
}

extern "C" {

void waxpby(int n, double *w, double alpha, const double *x, double beta, const double *y)
{
  waxpby_(n, w, alpha, x, beta, y);
}

double dot(int n, const double *x, const double *y)
{
  return dot_(n, x, y);
}

void pointwiseDivide(int n, double *w, const double *x, const double *y)
{
  pointwiseDivide_(n, w, x, y);
}

void pointwiseMultiply(int n, double *w, const double *x, const double *y)
{
  pointwiseMultiply_(n, w, x, y);
}

void waxpb(int n, double *w, double alpha, const double *x, double beta)
{
  waxpb_(n, w, alpha, x, beta);
}

double sum(int n, const double *x)
{
  double sum = 0;
#pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; i++) {
    sum += x[i];
  }
  return sum;
}

double minimum(int n, const double *x)
{
  double minimum = DBL_MAX;
#pragma omp parallel for reduction(min:minimum)
  for (int i = 0; i < n; i++) {
    minimum = min(minimum, x[i]);
  }
  return minimum;
}

void min(int n, double *w, const double *x, double alpha)
{
#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    w[i] = min(x[i], alpha);
  }
}

void CSR_abs(int n, double *w, const double *x)
{
#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    w[i] = fabs(x[i]);
  }
}

void CSR_exp(int n, double *w, const double *x)
{
#pragma omp parallel
  {
    int begin, end;
    SpMP::getSimpleThreadPartition(&begin, &end, n);
    vdExp(end - begin, x + begin, w + begin);
  }
/*#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    w[i] = exp(x[i]);
  }*/
}

void CSR_log1p(int n, double *w, const double *x)
{
#pragma omp parallel
  {
    int begin, end;
    SpMP::getSimpleThreadPartition(&begin, &end, n);
    vdLog1p(end - begin, x + begin, w + begin);
  }
/*#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    w[i] = log1p(x[i]);
  }*/
}

} // extern "C"
