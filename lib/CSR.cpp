// In this file, we assume the input matrix is in CSR format, in which  
// rowPtr and colIdx are 0-based.

#include <cstdio>
#include <cassert>
#include <cstring>
#include <climits>
#include <vector>

#include <omp.h>

#include "CSR.hpp"

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
static Graph *constructBoostTaskGraph(const CSR& A)
{
  vector<pair<int, int> > edges;
  for (int i = 0; i < A.m; ++i) {
    for (int j = A.rowPtr[i]; j < A.rowPtr[i + 1]; ++j) {
      int src = A.colIdx[j];
      if (src != i) {
        edges.push_back(make_pair(src, i));
      }
    }
  } // for each row i

  return new Graph(edges_are_unsorted_t(), edges.begin(), edges.end(), A.m);
}

void CSR::boostGetRCMPermutation(int *perm, int *inversePerm, int source /*=-1*/) const
{
  Graph *g = constructBoostTaskGraph(*this);

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

void CSR::permuteRowPtr_(CSR* out, const int *inversePerm) const
{
  out->rowPtr[0] = 0;

  int rowPtrSum[omp_get_max_threads() + 1];
  rowPtrSum[0] = 0;

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    int iPerThread = (m + nthreads - 1)/nthreads;
    int iBegin = min(iPerThread*tid, m);
    int iEnd = min(iBegin + iPerThread, m);

    out->rowPtr[iBegin] = 0;
    int i;
    for (i = iBegin; i < iEnd - 1; ++i) {
      int row = inversePerm ? inversePerm[i] : i;
      int begin = rowPtr[row], end = rowPtr[row + 1];
      out->rowPtr[i + 1] = out->rowPtr[i] + end - begin;
    }
    if (i < iEnd) {
      int row = inversePerm ? inversePerm[i] : i;
      int begin = rowPtr[row], end = rowPtr[row + 1];
      rowPtrSum[tid + 1] = out->rowPtr[i] + end - begin;
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
      out->rowPtr[m] = rowPtrSum[nthreads];
      assert(out->rowPtr[m] == rowPtr[m]);
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
    int begin = rowPtr[row], end = rowPtr[row + 1];
    int newBegin = out->rowPtr[i];

    int k = newBegin;
    for (int j = begin; j < end; ++j, ++k) {
      int c = colIdx[j];
      int newColIdx = columnPerm[c];

      out->colIdx[k] = newColIdx;
      out->values[k] = values[j];
    }

    memcpy(
      out->values + newBegin + (end - begin),
      values + end, (rowPtr[row + 1] - end)*sizeof(double));
    memcpy(
      out->colIdx + newBegin + (end - begin),
      colIdx + end, (rowPtr[row + 1] - end)*sizeof(int));

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
  int bw = INT_MIN;
  for (int i = 0; i < m; ++i) {
    for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
      int c = colIdx[j];
      int temp = c - i;
      if (temp < 0) temp = -temp;
      bw = max(temp, bw);
    }
  }
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
  for (int i = 0; i < rowPtr[m]; i++) {
    printf("%d ", colIdx[i]);
  }
  printf("\nValues: ");
  for (int i = 0; i < rowPtr[m]; i++) {
    printf("%f ", values[i]);
  }
  printf("\n\n");

  for (int i = 0; i < m; ++i) {
    int jj = 0;
    printf("%d: ", i);
    for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
      int c = colIdx[j];
      for ( ; jj < c; ++jj) printf("0 ");
      printf("%g ", values[j]);
      ++jj;
    }
    for ( ; jj < m; ++jj) printf("0 ");
    printf("\n");
  }
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

 void CSR::make0BasedIndexing() const
 {
#pragma omp parallel for
   for (int i = 0; i <= m; ++i) {
     rowPtr[i]--;
   }
#pragma omp parallel for
   for (int i = 0; i < rowPtr[m]; ++i) {
     colIdx[i]--;
   }
 }

 void CSR::make1BasedIndexing() const
 {
#pragma omp parallel for
   for (int i = 0; i < rowPtr[m]; ++i) {
     colIdx[i]++;
   }
#pragma omp parallel for
   for (int i = 0; i <= m; ++i) {
     rowPtr[i]++;
   }
 }
