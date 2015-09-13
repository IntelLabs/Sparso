#include <cstring>
#include <algorithm>

#include "IChol.hpp"
#include "SpMP/synk/barrier.hpp"

using namespace std;

namespace SpMP
{

void ichol0(CSR& A, double *l, const LevelSchedule& schedule)
{
  assert(A.isSymmetric());

  int base = A.getBase();

  const int *rowptr = A.rowptr - base;
  const int *colidx = A.colidx - base;
  const int *diagptr = A.diagptr - base;
  const double *values = A.values - base;

#ifndef NDEBUG
  for (int i = base; i < A.m + base; ++i) {
    assert(std::is_sorted(colidx + rowptr[i], colidx + rowptr[i + 1]));
    for (int j = rowptr[i]; j < diagptr[i]; ++j) {
      assert(colidx[j] < i);
    }
    assert(colidx[diagptr[i]] == i);
    for (int j = diagptr[i] + 1; j < rowptr[i + 1]; ++j) {
      assert(colidx[j] > i);
    }
  }
#endif

  l -= base;

  for (int i = base; i < A.getNnz() + base; ++i) {
    l[i] = values[i];
  }

  for (int i = base; i < A.m + base; ++i) {
    double a_ii = l[diagptr[i]] = sqrt(l[diagptr[i]]);

    for (int j = diagptr[i] + 1; j < rowptr[i + 1]; ++j) {
      int c = colidx[j];
      double tmp = l[j] /= a_ii;

      int k1 = j, k2 = diagptr[c];
      while (k1 < rowptr[i + 1] && k2 < rowptr[c + 1]) {
        if (colidx[k1] < colidx[k2]) ++k1;
        else if (colidx[k1] > colidx[k2]) ++k2;
        else {
          l[k2] -= tmp*l[k1];
          ++k1; ++k2;
        }
      }
    }
  } // for each row
}

} // namespace SpMP
