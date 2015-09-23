#include <cstring>

#include "SpMP/CSR.hpp"

#include "CSR_Interface.h"

using namespace std;

namespace SpMP
{

template<typename T, int BASE = 0>
void adb_(
  const T *A_values, const int *A_rowptr, const int *A_colidx, int A_m, int A_n,
  const T *B_values, const int *B_rowptr, const int *B_colidx, int B_m, int B_n,
  const double *d,
  T *C_data, const int *C_i, const int *C_j)
{
  if (A_n != B_m)
  {
    printf("Warning! incompatible matrix dimensions!\n");
    return;
  }

  T *C_temp_array = new T[B_n*omp_get_max_threads()];

#pragma omp parallel
  {
  int tid = omp_get_thread_num();
  T *C_temp = C_temp_array + tid*B_n;

#pragma omp for
  for (int i = 0; i < A_m; i++) {
    for (int j = C_i[i] - BASE; j < C_i[i + 1] - BASE; j++) {
      C_temp[C_j[j] - BASE] = 0;
    }
    for (int j = A_rowptr[i] - BASE; j < A_rowptr[i + 1] - BASE; j++) {
      int jcol = A_colidx[j] - BASE;
      T a_entry = A_values[j]*d[jcol];
      for (int k = B_rowptr[jcol] - BASE; k < B_rowptr[jcol + 1] - BASE; k++) {
        int kcol = B_colidx[k] - BASE;
        C_temp[kcol] += a_entry*B_values[k];
      }
    }
    for (int j = C_i[i] - BASE; j < C_i[i + 1] - BASE; j++) {
      C_data[j] = C_temp[C_j[j] - BASE];
    }
  } // for each row
  } // omp parallel

  delete[] C_temp_array;
} // adb_

template<typename T, int BASE = 0, bool SORT = true>
void inspectADB_(CSR * C, const CSR *A, const CSR *B)
{
  assert(A->getBase() == BASE);
  assert(B->getBase() == BASE);

  if (A->n != B->m)
  {
    printf("Warning! incompatible matrix dimensions!\n");
    return;
  }

  C->rowptr[0] = BASE;

  int *marker_array = new int[C->n*omp_get_max_threads()];

  int *marker = marker_array; 
  for (int i = 0; i < C->n; ++i) {
    marker[i] = -1;
  }

  int counter = 0;
  for (int i = 0; i < A->m; i++) {
    for (int j = A->rowptr[i] - BASE; j < A->rowptr[i + 1] - BASE; j++) {
      int jcol = A->colidx[j] - BASE;
      for (int k = B->rowptr[jcol] - BASE; k < B->rowptr[jcol + 1] - BASE; k++) {
        int kcol = B->colidx[k] - BASE;
        if (marker[kcol] != i) {
          marker[kcol] = i;
          ++counter;
        }
      }
    }
    C->rowptr[i + 1] = counter + BASE;
  } // for each row

  C->colidx = new int[counter];
  C->values = new T[counter];

  for (int i = 0; i < C->n; ++i) {
    marker[i] = -1;
  }

  counter = 0;
  for (int i = 0; i < A->m; i++) {
    for (int j = A->rowptr[i] - BASE; j < A->rowptr[i + 1] - BASE; j++) {
      int jcol = A->colidx[j] - BASE;
      for (int k = B->rowptr[jcol] - BASE; k < B->rowptr[jcol + 1] - BASE; k++) {
        int kcol = B->colidx[k] - BASE;
        if (marker[kcol] != i) {
          marker[kcol] = i;
          C->colidx[counter] = kcol + BASE;
          ++counter;
        }
      }
    }

    bool SORT = true;
    if (SORT) {
      for (int j = C->rowptr[i] + 1 - BASE; j < C->rowptr[i + 1] - BASE; j++) {
        int c = C->colidx[j];
        T v = C->values[j];

        int k = j - 1;
        while (k >= C->rowptr[i] - BASE && C->colidx[k] > c) {
          C->colidx[k + 1] = C->colidx[k];
          C->values[k + 1] = C->values[k];
          --k;
        }

        C->colidx[k + 1] = c;
        C->values[k + 1] = v;
      }
    }
  } // for each row

  delete[] marker;
}

void inspectADB(CSR *C, const CSR *A, const CSR *B)
{
  if (A->getBase() != B->getBase()) {
    fprintf(stderr, "Warning: matrix index bases don't match\n");
    return;
  }

  if (0 == A->getBase()) {
    inspectADB_<double, 0>(C, A, B);
  }
  else {
    assert(1 == A->getBase());
    inspectADB_<double, 1>(C, A, B);
  }
}

CSR *inspectADB(const CSR *A, const CSR *B)
{
  if (A->getBase() != B->getBase()) {
    fprintf(stderr, "Warning: matrix index bases don't match\n");
    return NULL;
  }

  CSR *C = new CSR;
  C->m = A->m;
  C->n = B->n;
  C->rowptr = new int[C->m + 1];

  inspectADB(C, A, B);

  return C;
}

void adb(CSR *C, const CSR *A, const CSR *B, const double *d)
{
  if (C->getBase() != A->getBase() || C->getBase() != B->getBase()) {
    fprintf(stderr, "Warning: matrix index bases don't match\n");
    return;
  }

  if (0 == C->getBase()) {
    adb_<double, 0>(
      A->values, A->rowptr, A->colidx, A->m, A->n,
      B->values, B->rowptr, B->colidx, B->m, B->n,
      d,
      C->values, C->rowptr, C->colidx);
  }
  else {
    assert(1 == C->getBase());
    adb_<double, 1>(
      A->values, A->rowptr, A->colidx, A->m, A->n,
      B->values, B->rowptr, B->colidx, B->m, B->n,
      d,
      C->values, C->rowptr, C->colidx);
  }
}
template<class T>
static void sort(int *col_idx, T *a, int start, int end)
{
  int i, j, it;
  double dt;

  for (i=end-1; i>start; i--)
    for(j=start; j<i; j++)
      if (col_idx[j] > col_idx[j+1]){

        if (a){
          dt=a[j]; 
          a[j]=a[j+1]; 
          a[j+1]=dt;
        }
        it=col_idx[j]; 
        col_idx[j]=col_idx[j+1]; 
        col_idx[j+1]=it;

      }
}

CSR *SpGEMMWithEps(CSR *A, CSR *B, double eps)
{
  double t = omp_get_wtime();
  double temp_t;

  int base = A->getBase();

  double *A_data = A->values - base;
  int *A_i = A->rowptr - base;
  int *A_j = A->colidx - base;
  int nrows_A = A->m;
  int ncols_A = A->n;

  double *B_data = B->values - base;
  int *B_i = B->rowptr - base;
  int *B_j = B->colidx - base;
  int nrows_B = B->m;
  int ncols_B = B->n;

  if (ncols_A != nrows_B)
  {
    printf("Warning! incompatible matrix dimensions!\n");
    return NULL;
  }

  int *C_i = MALLOC(int, nrows_A + 1);
  C_i -= base;

  int *C_j;
  double *C_data;

  int nthreads = omp_get_max_threads();
  int prefix_sum_workspace[nthreads + 1];
  unsigned long long flops = 0;

  int C_temp_size_per_thread = ((1 << 30) - nthreads - 1)/nthreads;
  int *C_temp_j_array = MALLOC(int, nthreads*C_temp_size_per_thread);
  double *C_temp_data_array = MALLOC(double, nthreads*C_temp_size_per_thread);
  double *x_array = MALLOC(double, nthreads*ncols_B);

#pragma omp parallel reduction(+:flops)
  {
    int ns, ne;
    SpMP::getSimpleThreadPartition(&ns, &ne, nrows_A);
    ns += base;
    ne += base;

    int ii = omp_get_thread_num();
    int num_threads = omp_get_num_threads();

    int *C_temp_j = C_temp_j_array + ii*C_temp_size_per_thread - base;
    double *C_temp_data = C_temp_data_array + ii*C_temp_size_per_thread - base;

    double *x = x_array + ii*ncols_B - base;
    for (int ib = base; ib < ncols_B + base; ib++)
      x[ib] = INFINITY;

    int num_nonzeros = base;
    for (int ic = ns; ic < ne; ic++)
    {
      int row_begin = C_i[ic] = num_nonzeros;
      for (int ia = A_i[ic]; ia < A_i[ic+1]; ia++)
      {
        int ja = A_j[ia];
        double a_entry = A_data[ia];
        for (int ib = B_i[ja]; ib < B_i[ja+1]; ib++)
        {
          int jb = B_j[ib];
          double b_entry = B_data[ib];
          if (INFINITY == x[jb])
          {
            x[jb] = a_entry*b_entry;
            C_temp_j[num_nonzeros] = jb;
            num_nonzeros++;
          }
          else
            x[jb] += a_entry*b_entry;
        }
      }
      int row_end = num_nonzeros;
      num_nonzeros = row_begin;
      for (int j = row_begin; j < row_end; ++j) {
        int jc = C_temp_j[j];
        if (jc == ic || fabs(x[jc]) > eps) {
          C_temp_data[num_nonzeros] = x[jc];
          C_temp_j[num_nonzeros] = jc;
          ++num_nonzeros;
        }
        x[jc] = INFINITY;
      }
      //sort(C_temp_j, C_temp_data, row_begin, num_nonzeros);
    }

    num_nonzeros -= base;

    int nnz;
    SpMP::prefixSum(&num_nonzeros, &nnz, prefix_sum_workspace);

#pragma omp master
    {
      C_i[nrows_A + base] = nnz + base;

      C_j = MALLOC(int, nnz);
      C_data = MALLOC(double, nnz);
    }
#pragma omp barrier

    for (int i = ns; i < ne; i++) {
      C_i[i] += num_nonzeros;
    }

    memcpy(C_j + num_nonzeros, C_temp_j + base, sizeof(int)*(prefix_sum_workspace[ii + 1] - prefix_sum_workspace[ii]));
    memcpy(C_data + num_nonzeros, C_temp_data + base, sizeof(double)*(prefix_sum_workspace[ii + 1] - prefix_sum_workspace[ii]));

    //if (0 == ii) {
      //printf("SpMP::SpGEMM phase 3 takes %f\n", omp_get_wtime() - temp_t);
    //}
  } /*end parallel region */

  FREE(C_temp_j_array);
  FREE(C_temp_data_array);
  FREE(x_array);

  C_i += base;

  //printf("SpMP::SpGEMM takes %f (nnz = %d, flops = %lld)\n", omp_get_wtime() - t, C_i[A->m + base] - base, flops);

  return new CSR(A->m, B->n, C_i, C_j, C_data);
}

CSR *SpAdd(double alpha, CSR *A, double beta, CSR *B)
{
  int base = A->getBase();

  double *A_data = A->values - base;
  int *A_i = A->rowptr - base;
  int *A_j = A->colidx - base;
  int nrows_A = A->m;
  int ncols_A = A->n;

  double *B_data = B->values - base;
  int *B_i = B->rowptr - base;
  int *B_j = B->colidx - base;
  int nrows_B = B->m;
  int ncols_B = B->n;

  if (nrows_A != nrows_B || ncols_A != ncols_B)
  {
    printf("Warning! incompatible matrix dimensions!\n");
    return NULL;
  }

  int *C_i = MALLOC(int, nrows_A + 1);
  C_i -= base;

  int *C_j;
  double *C_data;

  int nthreads = omp_get_max_threads();
  int prefix_sum_workspace[nthreads + 1];

  int C_temp_size_per_thread = ((1 << 30) - nthreads - 1)/nthreads;
  int *C_temp_j_array = MALLOC(int, nthreads*C_temp_size_per_thread);
  double *C_temp_data_array = MALLOC(double, nthreads*C_temp_size_per_thread);
  double *x_array = MALLOC(double, nthreads*ncols_B);

#pragma omp parallel
  {
    int ns, ne;
    SpMP::getSimpleThreadPartition(&ns, &ne, nrows_A);
    ns += base;
    ne += base;

    int ii = omp_get_thread_num();
    int num_threads = omp_get_num_threads();

    int *C_temp_j = C_temp_j_array + ii*C_temp_size_per_thread - base;
    double *C_temp_data = C_temp_data_array + ii*C_temp_size_per_thread - base;

    double *x = x_array + ii*ncols_B - base;
    for (int ib = base; ib < ncols_B + base; ib++)
      x[ib] = INFINITY;

    int num_nonzeros = base;
    for (int ic = ns; ic < ne; ic++)
    {
      int row_begin = C_i[ic] = num_nonzeros;
      for (int ia = A_i[ic]; ia < A_i[ic+1]; ia++)
      {
        int ja = A_j[ia];

        x[ja] = alpha*A_data[ia];
        C_temp_j[num_nonzeros] = ja;
        num_nonzeros++;
      }

      for (int ib = B_i[ic]; ib < B_i[ic+1]; ib++)
      {
        int jb = B_j[ib];
        double b_entry = beta*B_data[ib];
        if (INFINITY == x[jb])
        {
          x[jb] = b_entry;
          C_temp_j[num_nonzeros] = jb;
          num_nonzeros++;
        }
        else
          x[jb] += b_entry;
      }

      for (int j = row_begin; j < num_nonzeros; ++j) {
        int jc = C_temp_j[j];
        C_temp_data[j] = x[jc];
        x[jc] = INFINITY;
      }
      //sort(C_temp_j, C_temp_data, row_begin, num_nonzeros);
    }

    num_nonzeros -= base;

    int nnz;
    SpMP::prefixSum(&num_nonzeros, &nnz, prefix_sum_workspace);

#pragma omp master
    {
      C_i[nrows_A + base] = nnz + base;

      C_j = MALLOC(int, nnz);
      C_data = MALLOC(double, nnz);
    }
#pragma omp barrier

    for (int i = ns; i < ne; i++) {
      C_i[i] += num_nonzeros;
    }

    memcpy(C_j + num_nonzeros, C_temp_j + base, sizeof(int)*(prefix_sum_workspace[ii + 1] - prefix_sum_workspace[ii]));
    memcpy(C_data + num_nonzeros, C_temp_data + base, sizeof(double)*(prefix_sum_workspace[ii + 1] - prefix_sum_workspace[ii]));
  } /*end parallel region */

  FREE(C_temp_j_array);
  FREE(C_temp_data_array);
  FREE(x_array);

  C_i += base;

  return new CSR(A->m, B->n, C_i, C_j, C_data);
}

} // namespace SpMP
