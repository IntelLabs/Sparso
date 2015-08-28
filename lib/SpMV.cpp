#include <algorithm>

#include "SpMP/CSR.hpp"

#include "CSR_Interface.h"

using namespace std;
using namespace SpMP;

extern "C" {

void CSR_MultiplyWithVector(
  double *w,
  double alpha, const CSR_Handle *A, const double *x,
  double beta, const double *y,
  double gamma)
{
  ((CSR *)A)->multiplyWithVector(w, alpha, x, beta, y, gamma);
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
