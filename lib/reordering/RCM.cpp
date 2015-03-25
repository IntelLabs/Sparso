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
  memset(local_count, 0, sizeof(int)*omp_get_max_threads()*numLevels);

/*#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();*/

//#pragma omp for
    for (int i = 0; i < A->m; ++i) {
      ++local_count[/*tid*numLevels + */levels[i]];
    }

/*#pragma omp for
    for (int i = 0; i < levels; ++i) {
      for (int j = 1; j < nthreads; ++j) {
        local_count[i] += local_count[tid*numLevels + i];
      }
    }*/

//#pragma omp master
    {
      prefixSum[0] = 0;
      for (int i = 1; i <= numLevels; ++i) {
        prefixSum[i] = prefixSum[i - 1] + local_count[i - 1];
      }
    }
  //}

  delete[] local_count;
  assert(prefixSum[numLevels] == A->m);
}

int bfs(const CSR *A, int source, int *levels) {
  int numLevels = 0;
  levels[source] = numLevels;
  numLevels++;

  vector<int> q, nextQ;
  q.push_back(source);

  while (!q.empty()) {
    for (int u : q) {
      int index = nextQ.size();
      assert(levels[u] == numLevels - 1);

      for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
        int v = A->colIdx[j];
        if (levels[v] > levels[u] + 1) {
          levels[v] = levels[u] + 1;
          nextQ.push_back(v);
        }
      }
    } // for each u in current level

    q.clear();
    q.swap(nextQ);
    ++numLevels;
  } // while !q.empty()

  return numLevels;
}

void CSR::getRCMPermutationNew(int *perm, int *inversePerm, int source /*=-1*/) const
{
  // 1. Start vertex
  // TODO: pseudo diameter heuristic
  assert(source >= 0 && source < m);

  // 2. BFS
  int *levels = new int[m];
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
