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
  numLevels++;

  int *q = new int[A->m];
  int *nextQ = new int[A->m];
  q[0] = source;
  int qTail = 1;
  int nextQTail;

  while (qTail) {
    nextQTail = 0;
    for (int i = 0; i < qTail; ++i) {
      int u = q[i];
      assert(levels[u] == numLevels - 1);

      for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
        int v = A->colIdx[j];
        if (levels[v] == INT_MAX) {
          levels[v] = levels[u] + 1;
          nextQ[nextQTail] = v;
          ++nextQTail;
        }
      }
    } // for each u in current level

    swap(q, nextQ);
    qTail = nextQTail;
    ++numLevels;
  } // while q not empty

  delete[] q;
  delete[] nextQ;

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
