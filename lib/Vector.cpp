#include <cstdio>
#include <cassert>
#include <algorithm>

#include <omp.h>

#include "SpMP/Utils.hpp"

using namespace std;

// w = alpha*x + beta*w
template<class T>
void waxpby_(int n, T *w, T alpha, const T *x, T beta, const T *y)
{
  if (1 == alpha) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = x[i] + beta*y[i];
    }
  }
  else if (1 == beta) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i] + y[i];
    }
  }
  else {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = alpha*x[i] + beta*y[i];
    }
  }
}

template<class T>
void waxpb_(int n, T *w, T alpha, const T *x, T beta)
{
  if (1 == alpha) {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      w[i] = x[i] + beta;
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

void reorderVector(double *v, double *tmp, const int *perm, int len)
{
  return SpMP::reorderVector(v, tmp, perm, len);
}

void reorderVectorWithInversePerm(double *v, double *tmp, const int *inversePerm, int len)
{
  return SpMP::reorderVectorWithInversePerm(v, tmp, inversePerm, len);
}

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

} // extern "C"
