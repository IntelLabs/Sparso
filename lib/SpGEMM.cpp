#include "SpMP/CSR.hpp"

#include "CSR_Interface.h"

using namespace std;
using namespace SpMP;

template<typename T>
void adb_(
  const T *A_data, const int *A_i, const int *A_j, int nrows_A, int ncols_A,
  const T *B_data, const int *B_i, const int *B_j, int nrows_B, int ncols_B,
  const double *d,
  T *C_data, int *C_i, int *C_j, int *cnnz)
{
  int         ia, ib, ic, ja, jb, num_nonzeros=0;
  int	       row_start, counter;
  T              a_entry, b_entry;
  int         *B_marker;

  if (ncols_A != nrows_B)
  {
    printf("Warning! incompatible matrix dimensions!\n");
    return;
  }

  B_marker = new int[ncols_B];

  for (ib = 0; ib < ncols_B; ib++)
    B_marker[ib] = -1;

  C_i[0] = 0;
  for (ic = 0; ic < nrows_A; ic++)
  {
    for (ia = A_i[ic]; ia < A_i[ic+1]; ia++)
    {
      ja = A_j[ia];
      for (ib = B_i[ja]; ib < B_i[ja+1]; ib++)
      {
        jb = B_j[ib];
        if (B_marker[jb] != ic)
        {
          B_marker[jb] = ic;
          num_nonzeros++;
        }
      }
    }
    C_i[ic+1] = num_nonzeros;
  }

  *cnnz = num_nonzeros;

  for (ib = 0; ib < ncols_B; ib++)
    B_marker[ib] = -1;

  counter = 0;
  for (ic = 0; ic < nrows_A; ic++)
  {
    row_start = C_i[ic];
    for (ia = A_i[ic]; ia < A_i[ic+1]; ia++)
    {
      ja = A_j[ia];
      a_entry = A_data[ia]*d[ja];
      for (ib = B_i[ja]; ib < B_i[ja+1]; ib++)
      {
        jb = B_j[ib];
        b_entry = B_data[ib];
        if (B_marker[jb] < row_start)
        {
          B_marker[jb] = counter;
          C_j[B_marker[jb]] = jb;
          C_data[B_marker[jb]] = a_entry*b_entry;
          counter++;
        }
        else {
          C_data[B_marker[jb]] += a_entry*b_entry;
        }

      }
    }

    bool SORT = true;
    if (SORT) {
      for (int j = row_start + 1; j < counter; ++j) {
        int c = C_j[j];
        double v = C_data[j];

        int k = j - 1;
        while (k >= row_start && C_j[k] > c) {
          C_j[k + 1] = C_j[k];
          C_data[k + 1] = C_data[k];
          --k;
        }

        C_j[k + 1] = c;
        C_data[k + 1] = v;
      }
    } // SORT
  } // for each row

  delete[] B_marker;
} // adb_

extern "C" {

void CSR_ADB(CSR_Handle *C, const CSR_Handle *A, const CSR_Handle *B, const double *d)
{
  CSR *Ccsr = (CSR *)C;
  int nnzEstimated = Ccsr->rowptr[Ccsr->m] - 1;

  CSR *Acsr = (CSR *)A;
  CSR *Bcsr = (CSR *)B;

  Acsr->make0BasedIndexing();
  Bcsr->make0BasedIndexing();
  Ccsr->make0BasedIndexing();

  int nnzActual = -1;
  adb_(
    Acsr->values, Acsr->rowptr, Acsr->colidx, Acsr->m, Acsr->n,
    Bcsr->values, Bcsr->rowptr, Bcsr->colidx, Bcsr->m, Bcsr->n,
    d,
    Ccsr->values, Ccsr->rowptr, Ccsr->colidx, &nnzActual);

  Acsr->make1BasedIndexing();
  Bcsr->make1BasedIndexing();
  Ccsr->make1BasedIndexing();

  assert(nnzActual == nnzEstimated);
}

}
