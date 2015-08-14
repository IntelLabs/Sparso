#include "SpMP/CSR.hpp"

#include "CSR_Interface.h"

using namespace std;
using namespace SpMP;

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
  for (int i = 0; i < A_m; i++)
  {
    for (int j = C_i[i] - BASE; j < C_i[i + 1] - BASE; j++) {
      C_temp[C_j[j] - BASE] = 0;
    }
    for (int j = A_rowptr[i] - BASE; j < A_rowptr[i + 1] - BASE; j++)
    {
      int jcol = A_colidx[j] - BASE;
      T a_entry = A_values[j]*d[jcol];
      for (int k = B_rowptr[jcol] - BASE; k < B_rowptr[jcol + 1] - BASE; k++)
      {
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

extern "C" {

CSR_Handle *CSR_ADBInspect(const CSR_Handle *A, const CSR_Handle *B)
{
  const int BASE = 1;

  CSR *Acsr = (CSR *)A;
  CSR *Bcsr = (CSR *)B;

  assert(Acsr->base == BASE);
  assert(Bcsr->base == BASE);

  if (Acsr->n != Bcsr->m)
  {
    printf("Warning! incompatible matrix dimensions!\n");
    return NULL;
  }

  CSR *C = new CSR;
  C->m = Acsr->m;
  C->n = Bcsr->n;
  C->rowptr = new int[C->m + 1];
  C->rowptr[0] = BASE;
  C->base = BASE;

  int *marker_array = new int[C->n*omp_get_max_threads()];

  int *marker = marker_array; 
  for (int i = 0; i < C->n; ++i) {
    marker[i] = -1;
  }

  int counter = 0;
  for (int i = 0; i < Acsr->m; i++) {
    for (int j = Acsr->rowptr[i] - BASE; j < Acsr->rowptr[i + 1] - BASE; j++) {
      int jcol = Acsr->colidx[j] - BASE;
      for (int k = Bcsr->rowptr[jcol] - BASE; k < Bcsr->rowptr[jcol + 1] - BASE; k++) {
        int kcol = Bcsr->colidx[k] - BASE;
        if (marker[kcol] != i) {
          marker[kcol] = i;
          ++counter;
        }
      }
    }
    C->rowptr[i + 1] = counter + BASE;
  } // for each row

  C->colidx = new int[counter];
  C->values = new double[counter];

  for (int i = 0; i < C->n; ++i) {
    marker[i] = -1;
  }

  counter = 0;
  for (int i = 0; i < Acsr->m; i++) {
    for (int j = Acsr->rowptr[i] - BASE; j < Acsr->rowptr[i + 1] - BASE; j++) {
      int jcol = Acsr->colidx[j] - BASE;
      for (int k = Bcsr->rowptr[jcol] - BASE; k < Bcsr->rowptr[jcol + 1] - BASE; k++) {
        int kcol = Bcsr->colidx[k] - BASE;
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
        double v = C->values[j];

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

  return (CSR_Handle *)C;
}

void CSR_ADB(CSR_Handle *C, const CSR_Handle *A, const CSR_Handle *B, const double *d)
{
  CSR *Ccsr = (CSR *)C;
  CSR *Acsr = (CSR *)A;
  CSR *Bcsr = (CSR *)B;

  assert(Ccsr->base == 1);
  assert(Acsr->base == 1);
  assert(Bcsr->base == 1);

  adb_<double, 1>(
    Acsr->values, Acsr->rowptr, Acsr->colidx, Acsr->m, Acsr->n,
    Bcsr->values, Bcsr->rowptr, Bcsr->colidx, Bcsr->m, Bcsr->n,
    d,
    Ccsr->values, Ccsr->rowptr, Ccsr->colidx);
}

} // extern "C"
