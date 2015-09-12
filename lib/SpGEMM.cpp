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

} // namespace SpMP
