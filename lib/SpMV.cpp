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

}
