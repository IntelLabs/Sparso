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

  //void transpose(CSR *out) const;

  /**
   * get reverse Cuthill Mckee permutation that tends to reduce the bandwidth
   *
   * @param source starting vertex (-1 to use pseudo diameter heuristic)
   */
  void boostGetRCMPermutation(int *perm, int *inversePerm, int source = -1) const;
  void getRCMPermutation(int *perm, int *inversePerm, int source = -1) const;

  void permute(CSR *out, const int *columnPerm, const int *rowInversePerm) const;

  int getBandwidth() const;

  void make0BasedIndexing() const; // assume it's originally in 1-based indexing. This function is not idempotent
  void make1BasedIndexing() const; // assume it's originally in 0-based indexing.

  void printInDense() const;
  void printSomeValues(int distance) const;

  int m, n;
  int *rowPtr; // rowptr. i[0] = 0, i[m] = nnz
  int *colIdx; // colidx
  double *values;

protected:
  void permuteRowPtr_(CSR *out, const int *inversePerm) const;
}; // class CSR

#endif // CSR_HPP
