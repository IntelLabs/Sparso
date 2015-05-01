// In this file, we assume the input matrix is in CSR format, in which  
// rowPtr and colIdx are 0-based.

#include <cstdio>
#include <cassert>
#include <cstring>
#include <climits>
#include <vector>
#include <algorithm>

#include <omp.h>

#include "CSR.hpp"

#ifdef SEP
#include <sampling.h>
#endif

using namespace std;

void getInversePerm(int *inversePerm, const int *perm, int n)
{
#pragma omp parallel for
#pragma simd
  for (int i = 0; i < n; ++i) {
    inversePerm[perm[i]] = i;
  }
}

#ifdef USE_BOOST
#define BOOST_GRAPH_USE_NEW_CSR_INTERFACE
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/bandwidth.hpp>

using namespace boost;

typedef compressed_sparse_row_graph<directedS> Graph;

// convert CSR to boost Graph
template<int BASE = 0>
static Graph *constructBoostTaskGraph(const CSR& A)
{
  vector<pair<int, int> > edges;
  for (int i = 0; i < A.m; ++i) {
    for (int j = A.rowPtr[i] - BASE; j < A.rowPtr[i + 1] - BASE; ++j) {
      int src = A.colIdx[j] - BASE;
      if (src != i) {
        edges.push_back(make_pair(src, i));
      }
    }
  } // for each row i

  return new Graph(edges_are_unsorted_t(), edges.begin(), edges.end(), A.m);
}

void CSR::boostGetRCMPermutation(int *perm, int *inversePerm, int source /*=-1*/) const
{
  Graph *g = NULL;
  if (0 == base) {
    g = constructBoostTaskGraph(*this);
  }
  else {
    g = constructBoostTaskGraph<1>(*this);
  }

  vector<int> v;
  if (-1 == source) {
    cuthill_mckee_ordering(*g, back_inserter(v));
  }
  else {
    typedef out_degree_property_map<Graph> DegreeMap;
    std::vector<default_color_type> colors(num_vertices(*g));

    cuthill_mckee_ordering(*g, source, back_inserter(v), make_iterator_property_map(&colors[0], get(vertex_index, *g), colors[0]), make_out_degree_map(*g));
  }
  delete g;

  reverse_copy(v.begin(), v.end(), inversePerm);

  getInversePerm(perm, inversePerm, m);
}
#endif // USE_BOOST

void permuteRowPtr_(CSR* out, const CSR *in, const int *inversePerm)
{
  int m = in->m;

  out->base = in->base;
  out->rowPtr[0] = in->base;

  int rowPtrSum[omp_get_max_threads() + 1];
  rowPtrSum[0] = 0;

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    int iPerThread = (m + nthreads - 1)/nthreads;
    int iBegin = min(iPerThread*tid, m);
    int iEnd = min(iBegin + iPerThread, m);

    out->rowPtr[iBegin] = in->base;
    int i;
    for (i = iBegin; i < iEnd - 1; ++i) {
      int row = inversePerm ? inversePerm[i] : i;
      int begin = in->rowPtr[row], end = in->rowPtr[row + 1];
      out->rowPtr[i + 1] = out->rowPtr[i] + end - begin;
    }
    if (i < iEnd) {
      int row = inversePerm ? inversePerm[i] : i;
      int begin = in->rowPtr[row], end = in->rowPtr[row + 1];
      rowPtrSum[tid + 1] = out->rowPtr[i] + end - begin - in->base;
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
      out->rowPtr[m] = rowPtrSum[nthreads] + in->base;
      assert(out->rowPtr[m] == in->rowPtr[m]);
    }

    for (i = iBegin; i < iEnd; ++i) {
      out->rowPtr[i] += rowPtrSum[tid];
    }
  } // omp parallel
}

  /**
   * Compute y = A*x
   */
// TODO: remove this once MKL libray call is fine, or when reusing 
// works so that we can convert 0 to 1 based only once in the loop
// This is a temporary workaround. To remove in future.
template<int BASE = 0>
void CSR_MultiplyWithVector(int num_rows, const int *rowPtr, const int *colIdx, const double* values, const double *x, double *y)
{
//#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
#endif

#ifdef SEP
  static int cnt = 0;
  if (2 == cnt) {
    VTResumeSampling();
  }
#endif

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    int nnz = rowPtr[num_rows] - BASE;
    int nnzPerThread = (nnz + nthreads - 1)/nthreads;
    int iBegin = lower_bound(rowPtr, rowPtr + num_rows, nnzPerThread*tid + BASE) - rowPtr;
    int iEnd = lower_bound(rowPtr, rowPtr + num_rows, nnzPerThread*(tid + 1) + BASE) - rowPtr;
    assert(iBegin <= iEnd);
    assert(iBegin >= 0 && iBegin <= num_rows);
    assert(iEnd >= 0 && iEnd <= num_rows);

    for (int i = iBegin; i < iEnd; ++i) {
      double sum = 0;
      for (int j = rowPtr[i] - BASE; j < rowPtr[i + 1] - BASE; ++j) {
        sum += values[j]*x[colIdx[j] - BASE];
      }
      y[i] = sum;
    }
#ifdef MEASURE_LOAD_BALANCE
    double t = omp_get_wtime();
#pragma omp barrier
    barrierTimes[tid] = omp_get_wtime() - t;

#pragma omp barrier
#pragma omp master
    {
      double tEnd = omp_get_wtime();
      double barrierTimeSum = 0;
      for (int i = 0; i < nthreads; ++i) {
        barrierTimeSum += barrierTimes[i];
      }
      printf("%f load imbalance = %f\n", tEnd - tBegin, barrierTimeSum/(tEnd - tBegin)/nthreads);
    }
#undef MEASURE_LOAD_BALANCE
#endif // MEASURE_LOAD_BALANCE
  } // omp parallel

#ifdef SEP
  if (2 == cnt) {
    VTPauseSampling();
  }
  ++cnt;
#endif
}

extern "C" void CSR_MultiplyWithVector_1Based(int num_rows, int *rowPtr, int *colIdx, double* values, double *x, double *y)
{
  return CSR_MultiplyWithVector<1>(num_rows, rowPtr, colIdx, values, x, y);
}

void CSR::multiplyWithVector(double *y, const double *x) const
{
  return CSR_MultiplyWithVector<0>(m, rowPtr, colIdx, values, x, y);
}

template<int BASE = 0, bool SORT = false>
void permute_(
  CSR *out, const CSR *in,
  const int *columnPerm, const int *rowInversePerm)
{
  assert(in->base == BASE);
  permuteRowPtr_(out, in, rowInversePerm);

  int m = in->m;

//#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
#endif

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    int nnz = in->rowPtr[m] - BASE;
    int nnzPerThread = (nnz + nthreads - 1)/nthreads;
    int iBegin = lower_bound(out->rowPtr, out->rowPtr + m, nnzPerThread*tid + BASE) - out->rowPtr;
    int iEnd = lower_bound(out->rowPtr, out->rowPtr + m, nnzPerThread*(tid + 1) + BASE) - out->rowPtr;
    assert(iBegin <= iEnd);
    assert(iBegin >= 0 && iBegin <= m);
    assert(iEnd >= 0 && iEnd <= m);

    for (int i = iBegin; i < iEnd; ++i) {
      int row = rowInversePerm ? rowInversePerm[i] : i;
      int begin = in->rowPtr[row] - BASE, end = in->rowPtr[row + 1] - BASE;
      int newBegin = out->rowPtr[i] - BASE;

      int k = newBegin;
      for (int j = begin; j < end; ++j, ++k) {
        int c = in->colIdx[j] - BASE;
        int newColIdx = columnPerm[c];

        out->colIdx[k] = newColIdx + BASE;
        out->values[k] = in->values[j];
      }

      if (SORT) {
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
      }
    } // for each row

#ifdef MEASURE_LOAD_BALANCE
    double t = omp_get_wtime();
#pragma omp barrier
    barrierTimes[tid] = omp_get_wtime() - t;

#pragma omp barrier
#pragma omp master
    {
      double tEnd = omp_get_wtime();
      double barrierTimeSum = 0;
      for (int i = 0; i < nthreads; ++i) {
        barrierTimeSum += barrierTimes[i];
        int iBegin = lower_bound(out->rowPtr, out->rowPtr + m, nnzPerThread*i + BASE) - out->rowPtr;
        int iEnd = lower_bound(out->rowPtr, out->rowPtr + m, nnzPerThread*(i + 1) + BASE) - out->rowPtr;
        printf("%d %d-%d %d %f\n", i, iBegin, iEnd, out->rowPtr[iEnd] - out->rowPtr[iBegin], barrierTimes[i]);
      }
      printf("%f load imbalance = %f\n", tEnd - tBegin, barrierTimeSum/(tEnd - tBegin)/nthreads);
    }
#undef MEASURE_LOAD_BALANCE
#endif // MEASURE_LOAD_BALANCE
  } // omp parallel
}

void CSR::permute(
  CSR *out, const int *columnPerm, const int *rowInversePerm, bool sort/*=false*/) const
{
  if (0 == base) {
    if (sort) {
      permute_<0, true>(out, this, columnPerm, rowInversePerm);
    }
    else {
      permute_<0, false>(out, this, columnPerm, rowInversePerm);
    }
  }
  else {
    assert(1 == base);
    if (sort) {
      permute_<1, true>(out, this, columnPerm, rowInversePerm);
    }
    else {
      permute_<1, false>(out, this, columnPerm, rowInversePerm);
    }
  }
}

template<int BASE = 0>
int getBandwidth_(const CSR *A)
{
  int bw = INT_MIN;
#pragma omp parallel for reduction(max:bw)
  for (int i = 0; i < A->m; ++i) {
    for (int j = A->rowPtr[i] - BASE; j < A->rowPtr[i + 1] - BASE; ++j) {
      int c = A->colIdx[j] - BASE;
      int temp = c - i;
      if (temp < 0) temp = -temp;
      bw = max(temp, bw);
    }
  }
  return bw;
}

int CSR::getBandwidth() const
{
  if (0 == base) {
    return getBandwidth_<0>(this);
  }
  else {
    assert(1 == base);
    return getBandwidth_<1>(this);
  }
}

void printInDense(const CSR *A)
{
  // Raw format
  printf("RowPtr: ");
  for (int i = 0; i <= A->m; i++) {
    printf("%d ", A->rowPtr[i]);
  }
  printf("\nColIdx: ");
  for (int i = 0; i < A->rowPtr[A->m] - base; i++) {
    printf("%d ", A->colIdx[i]);
  }
  printf("\nValues: ");
  for (int i = 0; i < A->rowPtr[A->m] - base; i++) {
    printf("%f ", A->values[i]);
  }
  printf("\n\n");

  for (int i = 0; i < A->m; ++i) {
    int jj = 0;
    printf("%d: ", i);
    for (int j = A->rowPtr[i] - base; j < A->rowPtr[i + 1] - base; ++j) {
      int c = A->colIdx[j] - base;
      for ( ; jj < c; ++jj) printf("0 ");
      printf("%g ", A->values[j]);
      ++jj;
    }
    for ( ; jj < A->m; ++jj) printf("0 ");
    printf("\n");
  }
}

// print the CSR matrix every distance elements, as well as the first and 
// the last 10 elements
void CSR::printSomeValues(int distance) const
{
  fflush(stdout);
  printf("CSR values:\n");
  int count = 0;
  int nnz = rowPtr[m] - base;
  for (int i = 0; i < m; ++i) {
    for (int j = rowPtr[i] - base; j < rowPtr[i + 1] - base; ++j) {
      int c = colIdx[j] - base;
      if ((count % distance) == 0 || count < 10 || nnz - count < 10) {
          printf("%d %d %g\n", i + 1, c + 1, values[j]);
          // always print in 1-based for easier compare with .mtx file
      }
      count++;
    }
  }
  fflush(stdout);
}

/**
 * idx = idx2*dim1 + idx1
 * -> ret = idx1*dim2 + idx2
 *        = (idx%dim1)*dim2 + idx/dim1
 */
static inline int transpose_idx(int idx, int dim1, int dim2)
{
  return idx%dim1*dim2 + idx/dim1;
}

#if 0
void CSR::transpose(CSR *out)
{
   const double *A_data = A->values;
   const int *A_i = A->rowPtr;
   const int *A_j = A->colIdx;
   int num_rowsA = A->m;
   int num_colsA = A->n;
   int num_nonzerosA = A->rowPtr[m];
   const double *data = A->values;

   int *AT_i = B->rowPtr;
   int *AT_j = B->colIdx;
   double *AT_data = B->values;

   /*--------------------------------------------------------------
    * First, ascertain that num_cols and num_nonzeros has been set.
    * If not, set them.
    *--------------------------------------------------------------*/

   if (! num_nonzerosA)
   {
      num_nonzerosA = A_i[num_rowsA];
   }

   if (num_rowsA && ! num_colsA)
   {
      HYPRE_Int max_col = -1;
      HYPRE_Int i, j;
      for (i = 0; i < num_rowsA; ++i)
      {
          for (j = A_i[i]; j < A_i[i+1]; j++)
          {
              if (A_j[j] > max_col)
                 max_col = A_j[j];
          }
      }
      num_colsA = max_col+1;
   }

   double t = MPI_Wtime();

   HYPRE_Int *bucket = hypre_CTAlloc(
    HYPRE_Int, num_colsA*omp_get_max_threads());

#define MIN(a, b) (((a) <= (b)) ? (a) : (b))

#pragma omp parallel
   {
   int nthreads = omp_get_num_threads();
   int tid = omp_get_thread_num();

   int iPerThread = (num_rowsA + nthreads - 1)/nthreads;
   int iBegin = MIN(iPerThread*tid, num_rowsA);
   int iEnd = MIN(iBegin + iPerThread, num_rowsA);

   HYPRE_Int i, j;
   for (i = 0; i < num_colsA; ++i) {
     bucket[tid*num_colsA + i] = 0;
   }

   // count the number of keys that will go into each bucket
   for (j = A_i[iBegin]; j < A_i[iEnd]; ++j) {
     int idx = A_j[j];
     bucket[tid*num_colsA + idx]++;
   }

   // prefix sum
 #pragma omp barrier

   for (i = tid*num_colsA + 1; i < (tid + 1)*num_colsA; ++i) {
     int transpose_i = transpose_idx(i, nthreads, num_colsA);
     int transpose_i_minus_1 = transpose_idx(i - 1, nthreads, num_colsA);

     bucket[transpose_i] += bucket[transpose_i_minus_1];
   }

 #pragma omp barrier
 #pragma omp single
   {
     for (i = 1; i < nthreads; ++i) {
       int j0 = num_colsA*i - 1, j1 = num_colsA*(i + 1) - 1;
       int transpose_j0 = transpose_idx(j0, nthreads, num_colsA);
       int transpose_j1 = transpose_idx(j1, nthreads, num_colsA);

       bucket[transpose_j1] += bucket[transpose_j0];
     }

     AT_i[num_colsA] = num_nonzerosA;
   }

   if (tid > 0) {
     int transpose_i0 = transpose_idx(num_colsA*tid - 1, nthreads, num_colsA);

     for (i = tid*num_colsA; i < (tid + 1)*num_colsA - 1; ++i) {
       int transpose_i = transpose_idx(i, nthreads, num_colsA);

       bucket[transpose_i] += bucket[transpose_i0];
     }
   }

 #pragma omp barrier

   if (data) {
      for (i = iEnd - 1; i >= iBegin; --i) {
        for (j = A_i[i + 1] - 1; j >= A_i[i]; --j) {
          int idx = A_j[j];
          --bucket[tid*num_colsA + idx];

          int offset = bucket[tid*num_colsA + idx];

          AT_data[offset] = A_data[j];
          AT_j[offset] = i;
        }
      }
   }
   else {
      for (i = iEnd - 1; i >= iBegin; --i) {
        for (j = A_i[i + 1] - 1; j >= A_i[i]; --j) {
          int idx = A_j[j];
          --bucket[tid*num_colsA + idx];

          int offset = bucket[tid*num_colsA + idx];

          AT_j[offset] = i;
        }
      }
   }

#pragma omp barrier

#pragma omp for
   for (i = 0; i < num_colsA; ++i) {
     AT_i[i] = bucket[i];
   }

   } // omp parallel

   hypre_TFree(bucket);

   return 0;
 }
#endif

void CSR::make0BasedIndexing()
{
  if (1 == base) {
#pragma omp parallel for
    for (int i = 0; i <= m; ++i) {
      rowPtr[i]--;
    }
#pragma omp parallel for
    for (int i = 0; i < rowPtr[m]; ++i) {
      colIdx[i]--;
    }

    base = 0;
  }
}

void CSR::make1BasedIndexing()
{
  if (0 == base) {
#pragma omp parallel for
    for (int i = 0; i < rowPtr[m]; ++i) {
      colIdx[i]++;
    }
#pragma omp parallel for
    for (int i = 0; i <= m; ++i) {
      rowPtr[i]++;
    }

    base = 1;
  }
}
