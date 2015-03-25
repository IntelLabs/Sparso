#include <cassert>
#include <climits>
#include <cstring>
#include <cstdio>

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
      assert(levels[i] != INT_MAX);
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

/**
 * @return -1 if shortcircuited
 */
int bfs(const CSR *A, int source, int *levels, int *width = NULL, int *shortCircuitWidth = NULL) {
#pragma omp parallel for
  for (int i = 0; i < A->m; ++i) {
    levels[i] = INT_MAX;
  }

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
      if (width) {
        *width = max(*width, qTailPrefixSum[nthreads]);
      }
      if (shortCircuitWidth) {
        if (qTailPrefixSum[nthreads] > *shortCircuitWidth) {
          numLevels = -1;
        }
      }
    }
#pragma omp barrier

    if (qTailPrefixSum[nthreads] == 0 || numLevels == -1) break;

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
            /*if (187905 == source && numLevels >= 1033) {
              printf("%d: %d->%d\n", numLevels, source, v);
            }*/
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

int getMinDegreeNode(const CSR *A)
{
  int local_min[omp_get_max_threads()];
  int local_min_idx[omp_get_max_threads()];
  int global_min_idx;

#pragma omp parallel
  {
  int nthreads = omp_get_num_threads();
  int tid = omp_get_thread_num();

  int temp_min = INT_MAX;
  int temp_min_idx = -1;

#pragma omp for
  for (int i = 0; i < A->m; ++i) {
    int degree = A->rowPtr[i + 1] - A->rowPtr[i];
    if (degree < temp_min) {
      temp_min = degree;
      temp_min_idx = i;
    }
  }

  local_min[tid] = temp_min;
  local_min_idx[tid] = temp_min_idx;

#pragma omp barrier
#pragma omp master
  {
    int global_min = INT_MAX;
    global_min_idx = -1;
    for (int i = 0; i < nthreads; ++i) {
      if (local_min[i] < global_min) {
        global_min = local_min[i];
        global_min_idx = local_min_idx[i];
      }
    }
    //printf("global_min = %d\n", global_min);
  }
  } // omp parallel

  return global_min_idx;
}

int getVerticesAtLevel(int *vertices, const int *levels, int level, int m)
{
  int idx = 0;
#pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    if (levels[i] == level) {
      int idx_save = __sync_fetch_and_add(&idx, 1);
      vertices[idx_save] = i;
    }
  }

  sort(vertices, vertices + idx);

  return idx;
}

int selectSourceWithPseudoDiameter(const CSR *A, int *levels)
{
  int s = getMinDegreeNode(A);
  printf("%d is the min degree node\n", s);
  int e = -1;

  int *candidates = new int[A->m];

  while (e == -1) {
    int diameter = bfs(A, s, levels);
    printf("diameter from %d is %d\n", s, diameter);
    int nCandidates = getVerticesAtLevel(candidates, levels, diameter - 2, A->m);

    // sort by vertex by ascending degree
    for (int i = 1; i < nCandidates; ++i) {
      int c = candidates[i];
      int degree = A->rowPtr[c + 1] - A->rowPtr[c];

      int j = i - 1;
      while (j >= 0 && A->rowPtr[candidates[j] + 1] - A->rowPtr[candidates[j]] > degree) {
        candidates[j + 1] = candidates[j];
        --j;
      }

      candidates[j + 1] = c;
    }

    // select first 5 that are not adjacent to any previously chosen vertex
    int outIdx = 1;
    for (int i = 1; i < nCandidates && outIdx < 5; ++i) {
      int u = candidates[i];
      bool adjacent = false;
      for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
        int v = A->colIdx[j];
        for (int k = 0; k < outIdx; ++k) {
          if (candidates[k] == v) {
            adjacent = true;
            break;
          }
        }
        if (adjacent) break;
      }
      if (!adjacent) {
        candidates[outIdx++] = u;
      }
    }

    printf("has %d candidates, brought down to %d\n", nCandidates, outIdx);

    nCandidates = outIdx;

    int minWidth = INT_MAX;
    for (int i = 0; i < nCandidates; ++i) {
      int width = INT_MIN;
      int newDiameter = bfs(A, candidates[i], levels, &width, &minWidth);
      printf("(diameter, width) from %d is (%d, %d)\n", candidates[i], newDiameter, width);
      if (-1 == newDiameter) { // short circuited
        continue;
      }
      else if (newDiameter > diameter && width < minWidth) {
        s = candidates[i];
        printf("select %d as the new starting point\n", s);
        e = -1;
        break;
      }
      else if (width < minWidth) {
        minWidth = width;
        e = candidates[i];
      }
    }
  }

  delete[] candidates;

  return s;
}

void CSR::getRCMPermutation(int *perm, int *inversePerm, int source /*=-1*/) const
{
  // 1. Start vertex
  // TODO: pseudo diameter heuristic
  double t;
  int *levels = new int[m];

  if (-1 == source) {
    t = omp_get_wtime();
    source = selectSourceWithPseudoDiameter(this, levels);
    printf("source selection takes %f\n", omp_get_wtime() - t);
    printf("source = %d\n", source);
  }
  assert(source >= 0 && source < m);

  // 2. BFS
  t = omp_get_wtime();
  int maxDegree = INT_MIN;
#pragma omp parallel for reduction(max:maxDegree)
  for (int i = 0; i < m; ++i) {
    maxDegree = max(maxDegree, rowPtr[i + 1] - rowPtr[i]);
  }

  int numLevels = bfs(this, source, levels);
  printf("numLevels = %d\n", numLevels);

  printf("bfs takes %f\n", omp_get_wtime() - t);
  t = omp_get_wtime();

  // 3. Reorder
  int *prefixSum = new int[numLevels + 1];
  prefixSumOfLevels(prefixSum, this, levels, numLevels);

  printf("prefix sum takes %f\n", omp_get_wtime() - t);
  t = omp_get_wtime();

  inversePerm[0] = source;

  volatile int *read_offset = new int[numLevels + 1];
  volatile int *write_offset = new int[numLevels + 1];
  read_offset[0] = 0;
  write_offset[0] = 1;

  int *children_array = new int[omp_get_max_threads()*maxDegree];

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    int *children = children_array + tid*maxDegree;

#pragma omp for
    for (int l = 1; l <= numLevels; ++l) {
      read_offset[l] = prefixSum[l];
      write_offset[l] = prefixSum[l];
    }

    for (int l = tid; l < numLevels; l += nthreads) {
      while (read_offset[l] != prefixSum[l + 1]) {
        while (read_offset[l] == write_offset[l]); // spin
        int u = inversePerm[read_offset[l]];
        ++read_offset[l];
        int childrenIdx = 0;
        for (int j = rowPtr[u]; j < rowPtr[u + 1]; ++j) {
          int v = colIdx[j];
          if (levels[v] == l + 1) {
            children[childrenIdx] = v;
            ++childrenIdx;
            levels[v] = -1;
          }
        }
        // sort increasing order of degree
        for (int i = 1; i < childrenIdx; ++i) {
          int c = children[i];
          int degree = rowPtr[c + 1] - rowPtr[c];

          int j = i - 1;
          while (j >= 0 && rowPtr[children[j] + 1] - rowPtr[children[j]] > degree) {
            children[j + 1] = children[j];
            --j;
          }

          children[j + 1] = c;
        }

        for (int i = 0; i < childrenIdx; ++i) {
          int c = children[i];
          inversePerm[write_offset[l + 1]] = c;
          ++write_offset[l + 1];
        }
      }
    } // for each level
  } // omp parallel
  delete[] read_offset;
  delete[] write_offset;
  delete[] children_array;
  printf("place takes %f\n", omp_get_wtime() - t);
  delete[] levels;
  delete[] prefixSum;

  int *temp = new int[m];
  reverse_copy(inversePerm, inversePerm + m, temp);
  copy(temp, temp + m, inversePerm);
  delete[] temp;

  getInversePerm(perm, inversePerm, n);
}
