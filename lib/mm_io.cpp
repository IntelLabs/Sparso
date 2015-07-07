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

using namespace std;
using namespace SpMP;

COO coo;

extern "C" {

// The following two functions are for Julia
void load_matrix_market_step1 (char *file, int *sizes, bool force_symmetric /*=false*/)
{
  load_matrix_market(file, coo, force_symmetric);
  sizes[0] = (int)coo.isSymmetric;
  sizes[1] = coo.m;
  sizes[2] = coo.n;
  sizes[3] = coo.nnz;
}

#define T double

void load_matrix_market_step2 (char *file, T *a, int *j, int *i, int *sizes, bool one_based_CSR)
{
  int m = sizes[1];
  int n = sizes[2];

  CSR csr(m, n, i, j, a, one_based_CSR ? 1 : 0);
  dcoo2crs(&coo, &csr, false /* don't create separate diag data*/);
}

void load_matrix_market (char *file, T **a, int **aj, int **ai, int *is_symmetric, int *am, int *an, int *annz, bool one_based_CSR /*=false*/, bool force_symmetric /*=false*/)
{
    // Read into a COO array (stored statically insded mm_read)
    int sizes[] = {0, 0, 0, 0};
    load_matrix_market_step1(file, sizes, force_symmetric);
    
    *is_symmetric = sizes[0];
    *am = sizes[1];
    *an = sizes[2];
    *annz = sizes[3];
    
    // Allocate space for the CSR array
    *a = (T *)_mm_malloc(sizeof(T)*(*annz), 64);
    *aj = (int *)_mm_malloc(sizeof(int)*(*annz), 64);
    *ai = (int *)_mm_malloc(sizeof(int)*(*am+1), 64);

    // Convert COO to CSR
    load_matrix_market_step2(file, *a, *aj, *ai, sizes, one_based_CSR);
}

}
