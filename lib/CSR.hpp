#ifndef CSR_HPP
#define CSR_HPP

/**
 * A CSR instance doesn't "own" data.
 * The user is responsible for deallocating memory.
 * Assume zero-based indexing.
 */
class CSR {
public :
  CSR(int numRows, int numCols, int *i, int *j, double *v) :
    m(numRows), n(numCols), rowPtr(i), colIdx(j), values(v)
  {
  }

  /**
   * Compute y = A*x
   */
  void multiplyWithVector(double *y, const double *x) const
  {
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
      double sum = 0;
      for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
        sum += values[j]*x[colIdx[j]];
      }
      y[i] = sum;
    }
  }

  /**
   * get reverse Cuthill Mckee permutation that tends to reduce the bandwidth
   */
  void getRCMPermutation(int *perm, int *inversePerm) const;

  void permute(CSR *out, const int *columnPerm, const int *rowInversePerm) const;

  int getBandwidth() const;

  void printInDense() const;

  int m, n;
  int *rowPtr; // rowptr. i[0] = 0, i[m] = nnz
  int *colIdx; // colidx
  double *values;

protected:
  void permuteRowPtr_(CSR *out, const int *inversePerm) const;
}; // class CSR

#endif // CSR_HPP
