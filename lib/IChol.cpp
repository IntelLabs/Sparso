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
