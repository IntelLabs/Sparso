// In this file, we assume the input matrix is in CSR format, in which  
// rowPtr and colIdx are 1-based. We convert their values to 0-based
// on the fly so that the underlying library keeps as 0-based.

#include <cstdio>
#include <vector>

#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/bandwidth.hpp>

#include <omp.h>

#include "CSR.hpp"

using namespace std;
using namespace boost;

typedef compressed_sparse_row_graph<directedS> Graph;

// convert CSR to boost Graph
static Graph *constructBoostTaskGraph(const CSR& A)
{
  vector<pair<int, int> > edges;
  for (int i = 0; i < A.m; ++i) {
    for (int j = A.rowPtr[i] - 1; j < A.rowPtr[i + 1] - 1; ++j) {
      int src = A.colIdx[j] - 1;
      if (src != i) {
        edges.push_back(make_pair(src, i));
      }
    }
  } // for each row i

  return new Graph(edges_are_unsorted_t(), edges.begin(), edges.end(), A.m);
}

static void getInversePerm(int *inversePerm, const int *perm, int n)
{
#pragma omp parallel for
#pragma simd
  for (int i = 0; i < n; ++i) {
    inversePerm[perm[i]] = i;
  }
}

void CSR::getRCMPermutation(int *perm, int *inversePerm) const
{
  Graph *g = constructBoostTaskGraph(*this);

  vector<int> v;
  cuthill_mckee_ordering(*g, back_inserter(v));
  delete g;

  reverse_copy(v.begin(), v.end(), inversePerm);

  getInversePerm(perm, inversePerm, m);
}

void CSR::permuteRowPtr_(CSR* out, const int *inversePerm) const
{
  out->rowPtr[0] = 1;

  int rowPtrSum[omp_get_max_threads() + 1];
  rowPtrSum[0] = 0;

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    int iPerThread = (m + nthreads - 1) / nthreads;
    int iBegin = min(iPerThread*tid, m);
    int iEnd = min(iBegin + iPerThread, m);

    out->rowPtr[iBegin] = 1;
    int i;
    for (i = iBegin; i < iEnd - 1; ++i) {
      int row = inversePerm ? inversePerm[i] : i;
      int begin = rowPtr[row], end = rowPtr[row + 1];
      out->rowPtr[i + 1] = out->rowPtr[i] + end - begin;
    }
    if (i < iEnd) {
      int row = inversePerm ? inversePerm[i] : i;
      int begin = rowPtr[row], end = rowPtr[row + 1];
      rowPtrSum[tid + 1] = out->rowPtr[i] + end - begin - 1;
    }
    else {
      rowPtrSum[tid + 1] = 0;
    }

#pragma omp barrier
#pragma omp single
    {
      for (int tid = 1; tid < nthreads; ++tid) {
        rowPtrSum[tid + 1] += rowPtrSum[tid];
      }
      out->rowPtr[m] = rowPtrSum[nthreads] + 1;
    }

    for (i = iBegin; i < iEnd; ++i) {
      out->rowPtr[i] += rowPtrSum[tid];
    }
  } // omp parallel
}

void CSR::permute(CSR *out, const int *columnPerm, const int *rowInversePerm) const
{
  permuteRowPtr_(out, rowInversePerm);

#pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    int row = rowInversePerm ? rowInversePerm[i] : i;
    int begin = rowPtr[row] - 1, end = rowPtr[row + 1] - 1;
    int newBegin = out->rowPtr[i] - 1;

    int k = newBegin;
    for (int j = begin; j < end; ++j, ++k) {
      int c = colIdx[j] - 1;
      int newColIdx = columnPerm[c];

      out->colIdx[k] = newColIdx + 1;
      out->values[k] = values[j];
    }

    memcpy(
      out->values + newBegin + (end - begin),
      values + end, (rowPtr[row + 1] -1 - end)*sizeof(double));
    memcpy(
      out->colIdx + newBegin + (end - begin),
      colIdx + end, (rowPtr[row + 1] -1 - end)*sizeof(int));

    // insertion sort
    for (int j = newBegin + 1; j < newBegin + (end - begin); ++j) {
      int c = out->colIdx[j];
      double v = out->values[j];

      int k = j - 1;
      while (k >= newBegin && out->colIdx[k] > c) {
        out->colIdx[k + 1] = out->colIdx[k];
        out->values[k + 1] = out->values[k];
        --k;
      }

      out->colIdx[k + 1] = c;
      out->values[k + 1] = v;
    }
  } // for each row
}

int CSR::getBandwidth() const
{
  Graph *g = constructBoostTaskGraph(*this);
  int bw = bandwidth(*g);
  delete g;
  return bw;
}

void CSR::printInDense() const
{
  // Raw format
  printf("RowPtr: ");
  for (int i = 0; i <= m; i++) {
    printf("%d ", rowPtr[i]);
  }
  printf("\nColIdx: ");
  for (int i = 0; i < rowPtr[m] - 1; i++) {
    printf("%d ", colIdx[i]);
  }
  printf("\nValues: ");
  for (int i = 0; i < rowPtr[m] - 1; i++) {
    printf("%f ", values[i]);
  }
  printf("\n\n");

  for (int i = 0; i < m; ++i) {
    int jj = 0;
    printf("%d: ", i);
    for (int j = rowPtr[i] - 1; j < rowPtr[i + 1] - 1; ++j) {
      int c = colIdx[j] - 1;
      for (; jj < c; ++jj) printf("0 ");
      printf("%g ", values[j]);
      ++jj;
    }
    for (; jj < m; ++jj) printf("0 ");
    printf("\n");
  }
}
