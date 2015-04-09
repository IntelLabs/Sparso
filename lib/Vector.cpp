extern "C" void reorderVector(double *v, double *tmp, const int *perm, int len)
{
  if (!perm) return;

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    tmp[perm[i]] = v[i];
  }

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    v[i] = tmp[i];
  }
}

extern "C" void reorderVectorWithInversePerm(double *v, double *tmp, const int *inversePerm, int len)
{
  if (!inversePerm) return;

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    tmp[i] = v[inversePerm[i]];
  }

#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    v[i] = tmp[i];
  }
}

  /**
   * Compute y = A*x
   */
// TODO: remove this once MKL libray call is fine, or when reusing 
// works so that we can convert 0 to 1 based only once in the loop
// This is a temporary workaround. To remove in future.
extern "C" void CSR_MultiplyWithVector_1Based(int num_rows, int *rowPtr, int *colIdx, double* values, double *x, double *y)
{
#pragma omp parallel for
    for (int i = 0; i < num_rows; ++i) {
      double sum = 0;
      for (int j = rowPtr[i] - 1; j < rowPtr[i + 1] - 1; ++j) {
        sum += values[j]*x[colIdx[j] - 1];
      }
      y[i] = sum;
    }
  }
