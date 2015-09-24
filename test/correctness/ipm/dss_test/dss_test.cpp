#include <string>
#include <mkl.h>

#include "SpMP/CSR.hpp"
#include "SpMP/Vector.hpp"

#include "CSR_Interface.h"

using namespace std;
using namespace SpMP;

// from mkl/solverc/source/dss_sym_c.c
#define NROWS       5
#define NCOLS       5
#define NNONZEROS   9
#define NRHS        1
static const MKL_INT nRows = NROWS;
static const MKL_INT nCols = NCOLS;
static const MKL_INT nNonZeros = NNONZEROS;
static const MKL_INT nRhs = NRHS;
static _INTEGER_t rowIndex[NROWS + 1] = { 1, 6, 7, 8, 9, 10 };
static _INTEGER_t columns[NNONZEROS] = { 1, 2, 3, 4, 5, 2, 3, 4, 5 };
static _DOUBLE_PRECISION_t values[NNONZEROS] = { 9, 1.5, 6, .75, 3, 0.5, 12, .625, 16 };

int main(int argc, const char *argv[])
{
  // Load matrix and vectors
  CSR A((string(argv[1]) + "-A.mtx").c_str());
  //CSR A((string(argv[1])).c_str());
  //printf("%d\n", A.m);
  //A.m = 1036;

  double *b, *p;
  int bLen, pLen;
  int dummy;
  loadVectorMatrixMarket((string(argv[1]) + "-b.mtx").c_str(), &b, &bLen, &dummy);
  assert(dummy == 1);
  loadVectorMatrixMarket((string(argv[1]) + "-p.mtx").c_str(), &p, &pLen, &dummy);
  assert(dummy == 1);

  // Compute A*A'
  CSR AT = *A.transpose();
  A.make1BasedIndexing();
  AT.make1BasedIndexing();
  CSR ADAT = *((CSR *)CSR_ADBInspect((CSR_Handle *)&A, (CSR_Handle *)&AT));

  double *d = MALLOC(double, A.n);
  fill(d, d + A.n, 1.);
  CSR_ADB((CSR_Handle *)&ADAT, (CSR_Handle *)&A, (CSR_Handle *)&AT, d);
  assert(ADAT.isSymmetric());
  //ADAT.print();

  //ADAT.m = nRows;
  //ADAT.n = nCols;
  //ADAT.rowptr = rowIndex;
  //ADAT.colidx = columns;
  //ADAT.values = values;

  _MKL_DSS_HANDLE_t handle;
  MKL_INT opt = MKL_DSS_MSG_LVL_WARNING;
  MKL_INT error = dss_create(handle, opt);
  if (error != MKL_DSS_SUCCESS) {
    printf("Solver returned error code %d\n", error);
    return -1;
  }

  MKL_INT sym = MKL_DSS_NON_SYMMETRIC;
  int nnz = ADAT.getNnz();
  error = dss_define_structure(
    handle, sym, ADAT.rowptr, ADAT.m, ADAT.n, ADAT.colidx, nnz);
  if (error != MKL_DSS_SUCCESS) {
    printf("Solver returned error code %d\n", error);
    return -1;
  }

  opt = MKL_DSS_AUTO_ORDER;
  error = dss_reorder(handle, opt, NULL);
  if (error != MKL_DSS_SUCCESS) {
    printf("dss_reorder returned error code %d\n", error);
    return -1;
  }

  MKL_INT type = MKL_DSS_POSITIVE_DEFINITE;
  error = dss_factor_real(handle, type, ADAT.values);
  if (error != MKL_DSS_SUCCESS) {
    printf("dss_factor_real returned error code %d\n", error);
    return -1;
  }

  double *x = MALLOC(double, ADAT.n);
  for (int i = 0; i < ADAT.n; i++) x[i] = i;
  mkl_dcsrgemv("n", &ADAT.m, ADAT.values, ADAT.rowptr, ADAT.colidx, x, b);
  for (int i = 0; i < ADAT.n; i++) x[i] = 0;

  opt = MKL_DSS_DEFAULTS;
  int nrhs = 1;
  error = dss_solve_real(handle, opt, b, nrhs, x);
  if (error != MKL_DSS_SUCCESS) {
    printf("dss_solve_real returned error code %d\n", error);
    return -1;
  }

  double eps = 1e-6;
  for (int i = 0; i < ADAT.n; i++) {
    if (x[i] > i + eps || x[i] < i - eps) {
      printf("Incorrect solution\n");
      return -1;
    }
  }

  return 0;
}
