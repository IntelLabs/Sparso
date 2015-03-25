#include <cassert>
#include <climits>
#include <cstring>

#include <vector>
#include <algorithm>

#include <omp.h>

#include "../CSR.hpp"

using namespace std;

void getInversePerm(int *inversePerm, const int *perm, int n);

// compute prefix sum of levels
static void prefixSumOfLevels(
  int *prefixSum,
  const CSR *A, const int *levels, int numLevels)
{
  int *local_count = new int[omp_get_max_threads()*numLevels];
  int *local_sum_array = new int[omp_get_max_threads() + 1];

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    memset(local_count + tid*numLevels, 0, sizeof(int)*numLevels);

#pragma omp for
    for (int i = 0; i < A->m; ++i) {
      ++local_count[tid*numLevels + levels[i]];
    }

    int lPerThread = (numLevels + nthreads - 1)/nthreads;
    int lBegin = min(lPerThread*tid, numLevels);
    int lEnd = min(lBegin + lPerThread, numLevels);

    int local_sum = 0;

    for (int l = lBegin; l < lEnd; ++l) {
      prefixSum[l] = local_sum;
      for (int t = 0; t < nthreads; ++t) {
        local_sum += local_count[t*numLevels + l];
      }
    }

    local_sum_array[tid + 1] = local_sum;

#pragma omp barrier
    if (0 == tid)
    {
      for (int t = 1; t <= nthreads; ++t) {
        local_sum_array[t + 1] += local_sum_array[t];
      }
      assert(local_sum_array[nthreads] == A->m);
      prefixSum[numLevels] = A->m;
    }
#pragma omp barrier

    if (tid > 0) {
      local_sum = local_sum_array[tid];
      for (int l = lBegin; l < lEnd; ++l) {
        prefixSum[l] += local_sum;
      }
    }
  } // omp parallel

  delete[] local_count;
  delete[] local_sum_array;
}

int bfs(const CSR *A, int source, int *levels) {
  int numLevels = 0;
  levels[source] = numLevels;

  int *q[2];
  q[0] = new int[A->m*omp_get_max_threads()];
  q[1] = new int[A->m*omp_get_max_threads()];
  q[0][0] = source;

  int *qTail = new int[omp_get_max_threads()];
  qTail[0] = 1;

  int *qTailPrefixSum = new int[omp_get_max_threads() + 1];
  qTailPrefixSum[0] = 0;

  int *rowPtrs = new int[omp_get_max_threads()*A->m];
  rowPtrs[0] = 0;
  rowPtrs[1] = A->rowPtr[source + 1] - A->rowPtr[source];

  int *nnzPrefixSum = new int[omp_get_max_threads() + 1];
  nnzPrefixSum[0] = 0;
  nnzPrefixSum[1] = rowPtrs[1];

#pragma omp parallel
  {
  int tid = omp_get_thread_num();
  int nthreads = omp_get_num_threads();

  if (tid > 0) {
    qTail[tid] = 0;
    nnzPrefixSum[tid + 1] = 0;
  }

  while (true) {
#pragma omp barrier
#pragma omp master
    {
      for (int t = 0; t < nthreads; ++t) {
        qTailPrefixSum[t + 1] = qTailPrefixSum[t] + qTail[t];
        nnzPrefixSum[t + 1] += nnzPrefixSum[t];
      }
      ++numLevels;
    }
#pragma omp barrier

    if (qTailPrefixSum[nthreads] == 0) break;

    // partition based on # of nnz
    int nnzPerThread = (nnzPrefixSum[nthreads] + nthreads - 1)/nthreads;
    int tBegin = upper_bound(
        nnzPrefixSum, nnzPrefixSum + nthreads + 1,
        nnzPerThread*tid) -
      nnzPrefixSum - 1;
    int tEnd = upper_bound(
        nnzPrefixSum, nnzPrefixSum + nthreads + 1,
        nnzPerThread*(tid + 1)) -
      nnzPrefixSum - 1;
    assert(tBegin >= 0 && tBegin <= nthreads);
    assert(tEnd >= 0 && tEnd <= nthreads);
    assert(tBegin <= tEnd);

    int iBegin, iEnd;
    if (0 == tid) {
      iBegin = 0;
    }
    else if (tBegin == nthreads) {
      iBegin = qTail[tEnd - 1];
    }
    else {
      iBegin = upper_bound(
          rowPtrs + tBegin*A->m, rowPtrs + tBegin*A->m + qTail[tBegin],
          nnzPerThread*tid - nnzPrefixSum[tBegin]) -
        (rowPtrs + tBegin*A->m) - 1;
    }

    if (tEnd == nthreads) {
      iEnd = 0;
    }
    else {
      iEnd = upper_bound(
          rowPtrs + tEnd*A->m, rowPtrs + tEnd*A->m + qTail[tEnd],
          nnzPerThread*(tid + 1) - nnzPrefixSum[tEnd]) -
        (rowPtrs + tEnd*A->m) - 1;
    }

#pragma omp barrier

    int *tailPtr = q[numLevels%2] + tid*A->m;
    int *rowPtr = rowPtrs + tid*A->m;
    *rowPtr = 0;

    for (int t = tBegin; t <= tEnd; ++t) {
      for (int i = (t == tBegin ? iBegin : 0);
          i < (t == tEnd ? iEnd : qTail[t]);
          ++i) {
        int u = q[1 - numLevels%2][t*A->m + i];
        assert(levels[u] == numLevels - 1);

        for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
          int v = A->colIdx[j];
          if (__sync_bool_compare_and_swap(levels + v, INT_MAX, numLevels)) {
            *tailPtr = v;
            *(rowPtr + 1) =
              *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

            ++tailPtr;
            ++rowPtr;
          }
        }
      } // for each current node u
    }

#pragma omp barrier

    qTail[tid] = tailPtr - (q[numLevels%2] + tid*A->m);
    nnzPrefixSum[tid + 1] = *rowPtr;
  } // while true
  } // omp parallel

  delete[] q[0];
  delete[] q[1];
  delete[] rowPtrs;
  delete[] qTail;
  delete[] qTailPrefixSum;
  delete[] nnzPrefixSum;

  return numLevels;
}

void CSR::getRCMPermutationNew(int *perm, int *inversePerm, int source /*=-1*/) const
{
  // 1. Start vertex
  // TODO: pseudo diameter heuristic
  assert(source >= 0 && source < m);

  // 2. BFS
  int *levels = new int[m];
#pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    levels[i] = INT_MAX;
  }

  int numLevels = bfs(this, source, levels);
  printf("numLevels = %d\n", numLevels);

  // 3. Reorder
  int *prefixSum = new int[numLevels + 1];
  prefixSumOfLevels(prefixSum, this, levels, numLevels);

  inversePerm[0] = source;
  int read_offset = 0;
  int write_offset = 1;

  for (int l = 0; l < numLevels; ++l) {
    while (read_offset != prefixSum[l + 1]) {
      int u = inversePerm[read_offset];
      ++read_offset;
      vector<int> children;
      for (int j = rowPtr[u]; j < rowPtr[u + 1]; ++j) {
        int v = colIdx[j];
        if (levels[v] == l + 1) {
          children.push_back(v);
          levels[v] = -1;
        }
      }
      sort(children.begin(), children.end());
      for (int c : children) {
        inversePerm[write_offset] = c;
        ++write_offset;
      }
    }
  } // for each level
  delete[] levels;

  int *temp = new int[m];
  reverse_copy(inversePerm, inversePerm + m, temp);
  copy(temp, temp + m, inversePerm);
  delete[] temp;

  getInversePerm(perm, inversePerm, n);
}
