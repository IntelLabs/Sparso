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
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <limits.h>
# include <time.h>
# include <assert.h>
#include <algorithm>
#include <omp.h>

#include "SpMP/COO.hpp"
#include "SpMP/mm_io.h"
#include "SpMP/Vector.hpp"

using namespace std;
using namespace SpMP;

COO coo;

extern "C" {

// The following two functions are for Julia
void load_matrix_market_step1(char *file, int *sizes, bool force_symmetric /*=false*/, bool transpose /*=false*/)
{
  if (transpose) {
    loadMatrixMarketTransposed(file, coo);
  }
  else {
    loadMatrixMarket(file, coo, force_symmetric);
  }
  sizes[0] = (int)coo.isSymmetric;
  sizes[1] = coo.m;
  sizes[2] = coo.n;
  sizes[3] = coo.nnz;
}

#define T double

void load_matrix_market_step2(
  char *file, int *rowptr, int *colidx, T *values, int *sizes, bool one_based_CSR)
{
  int m = sizes[1];
  int n = sizes[2];
  int nnz = sizes[3];

  dcoo2csr(m, nnz, rowptr, colidx, values, coo.rowidx, coo.colidx, coo.values);

  if (one_based_CSR) {
    CSR csr;
    csr.m = m;
    csr.n = n;
    csr.rowptr = rowptr;
    csr.colidx = colidx;
    csr.values = values;

    csr.make1BasedIndexing();
  }
}

void load_matrix_market(
  char *file,
  int **rowptr, int **colidx, T **values,
  int *is_symmetric, int *m, int *n, int *nnz,
  bool one_based_CSR /*=false*/, bool force_symmetric /*=false*/)
{
  // Read into a COO array (stored statically insded mm_read)
  int sizes[] = {0, 0, 0, 0};
  load_matrix_market_step1(file, sizes, force_symmetric, false /*no transpose*/);

  *is_symmetric = sizes[0];
  *m = sizes[1];
  *n = sizes[2];
  *nnz = sizes[3];

  // Allocate space for the CSR array
  *rowptr = (int *)_mm_malloc(sizeof(int)*(*m+1), 64);
  *colidx = (int *)_mm_malloc(sizeof(int)*(*nnz), 64);
  *values = (T *)_mm_malloc(sizeof(T)*(*nnz), 64);

  // Convert COO to CSR
  load_matrix_market_step2(file, *rowptr, *colidx, *values, sizes, one_based_CSR);
}

void load_vector_matrix_market(const char *fileName, double **v, int *m, int *n)
{
  SpMP::loadVectorMatrixMarket(fileName, v, m, n);
}

void store_matrix_market(char *file, int m, int n, int *colptr, int *rowval, double *nzval)
{
  CSR A(m, n, colptr, rowval, nzval);
  A.storeMatrixMarket(file);
}

} // extern "C"
