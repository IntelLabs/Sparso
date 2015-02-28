void reorderVector(double *v, double *tmp, const int *perm, int len)
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

void reorderVectorWithInversePerm(double *v, double *tmp, const int *inversePerm, int len)
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
