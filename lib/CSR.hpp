#ifndef CSR_HPP
#define CSR_HPP

/**
 * A CSR instance doesn't "own" data.
 * The user is responsible for deallocating memory.
 * Assume zero-based indexing.
 */
class CSR {
public :
  CSR(int numRows, int numCols, int *i, int *j, double *v, int base = 0) :
    m(numRows), n(numCols), rowPtr(i), colIdx(j), values(v), base(base)
  {
  }

  /**
   * Compute y = A*x
   */
  void multiplyWithVector(double *y, const double *x) const;

  //void transpose(CSR *out) const;

  /**
   * get reverse Cuthill Mckee permutation that tends to reduce the bandwidth
   *
   * @param source starting vertex (-1 to use pseudo diameter heuristic)
   */
  void getRCMPermutation(int *perm, int *inversePerm, int source = -1) const;

#ifdef USE_BOOST
  void boostGetRCMPermutation(int *perm, int *inversePerm, int source = -1) const;
#endif

  /**
   * @param sort true if we want to sort nnzs of each row based on colidx
   */
  void permute(CSR *out, const int *columnPerm, const int *rowInversePerm, bool sort = false) const;

  int getBandwidth() const;

  void make0BasedIndexing(); // assume it's originally in 1-based indexing. This function is not idempotent
  void make1BasedIndexing(); // assume it's originally in 0-based indexing.

  void printInDense() const;
  void printSomeValues(int distance, bool is_1_based) const;

  int m, n;
  int *rowPtr; // rowptr. rowPtr[0] = 0, rowPtr[m] = nnz
  int *colIdx; // colidx
  double *values;

  int base; // 0 : 0-based, 1-based

protected:
  void permuteRowPtr_(CSR *out, const int *inversePerm) const;
}; // class CSR

#endif // CSR_HPP
