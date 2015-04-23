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
  const CSR *A, const int *levels, int numLevels,
  const int *components, int numOfComponents)
{
  int *local_count = new int[omp_get_max_threads()*numLevels];
  int *local_sum_array = new int[omp_get_max_threads() + 1];

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    memset(local_count + tid*numLevels, 0, sizeof(int)*numLevels);

    if (NULL == components) {
#pragma omp for
      for (int i = 0; i < A->m; ++i) {
        assert(levels[i] != INT_MAX);
        ++local_count[tid*numLevels + levels[i]];
      }
    }
    else {
#pragma omp for
      for (int i = 0; i < numOfComponents; ++i) {
        assert(components[i] >= 0 && components[i] < A->m);
        assert(levels[components[i]] != INT_MAX);
        ++local_count[tid*numLevels + levels[components[i]]];
      }
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
      assert(local_sum_array[nthreads] == numOfComponents);
      prefixSum[numLevels] = numOfComponents;
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

struct bfsAuxData
{
  int *q[2];
  int *qTail;
  int *qTailPrefixSum;
  int *rowPtrs;
  int *nnzPrefixSum;
  int *candidates;

  bfsAuxData(int m);
  ~bfsAuxData();
};

bfsAuxData::bfsAuxData(int m)
{
  q[0] = new int[m*omp_get_max_threads()];
  q[1] = new int[m*omp_get_max_threads()];

  qTail = new int[omp_get_max_threads()];
  qTailPrefixSum = new int[omp_get_max_threads() + 1];

  rowPtrs = new int[omp_get_max_threads()*m];
  nnzPrefixSum = new int[omp_get_max_threads() + 1];

  candidates = new int[m];
}

bfsAuxData::~bfsAuxData()
{
  delete[] q[0];
  delete[] q[1];
  delete[] qTail;
  delete[] qTailPrefixSum;
  delete[] rowPtrs;
  delete[] nnzPrefixSum;
  delete[] candidates;
}

/**
 * @return -1 if shortcircuited num of levels otherwise
 *
 * pre-condition: levels should be initialized to -1
 */
template<bool OUTPUT_VISITED = false>
int bfs(
  const CSR *A, int source, int *levels,
  bfsAuxData *aux,
  int *visited = NULL, int *numOfVisited = NULL,
  int *width = NULL, int *shortCircuitWidth = NULL) {

  int numLevels = 0;
  levels[source] = numLevels;

  int **q = aux->q;
  q[0][0] = source;

  int *qTail = aux->qTail;
  qTail[0] = 1;

  int *qTailPrefixSum = aux->qTailPrefixSum;
  qTailPrefixSum[0] = 0;

  int *rowPtrs = aux->rowPtrs;
  rowPtrs[0] = 0;
  rowPtrs[1] = A->rowPtr[source + 1] - A->rowPtr[source];

  int *nnzPrefixSum = aux->nnzPrefixSum;
  nnzPrefixSum[0] = 0;
  nnzPrefixSum[1] = rowPtrs[1];

  if (OUTPUT_VISITED) *numOfVisited = 0;

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
      if (OUTPUT_VISITED) *numOfVisited += qTailPrefixSum[nthreads];
    }
#pragma omp barrier

    if (qTailPrefixSum[nthreads] == 0 || numLevels == -1) break;

    // partition based on # of nnz
    int nnzPerThread = (nnzPrefixSum[nthreads] + nthreads - 1)/nthreads;
    int tBegin = upper_bound(
        nnzPrefixSum, nnzPrefixSum + nthreads + 1,
        nnzPerThread*tid) -
      nnzPrefixSum - 1;
    if (0 == tid) {
      tBegin = 0;
    }
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

    if (OUTPUT_VISITED) {
      memcpy(
        visited + *numOfVisited - qTailPrefixSum[nthreads] + qTailPrefixSum[tid], 
        q[1 - numLevels%2] + tid*A->m,
        sizeof(int)*(qTailPrefixSum[tid + 1] - qTailPrefixSum[tid]));
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
          if (OUTPUT_VISITED) {
            if (__sync_bool_compare_and_swap(levels + v, INT_MAX, numLevels)) {
              *tailPtr = v;
              *(rowPtr + 1) = *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

              ++tailPtr;
              ++rowPtr;
              /*if (187905 == source && numLevels >= 1033) {
                printf("%d: %d->%d\n", numLevels, source, v);
              }*/
            }
          }
          else {
            if (INT_MAX == levels[v]) {
              levels[v] = numLevels;
              *tailPtr = v;
              *(rowPtr + 1) = *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

              ++tailPtr;
              ++rowPtr;
            }
          }
        }
      } // for each current node u
    }

#pragma omp barrier

    qTail[tid] = tailPtr - (q[numLevels%2] + tid*A->m);
    nnzPrefixSum[tid + 1] = *rowPtr;
  } // while true
  } // omp parallel

#ifndef NDEBUG
  if (OUTPUT_VISITED) {
    int *temp = new int[*numOfVisited];
    copy(visited, visited + *numOfVisited, temp);
    sort(temp, temp + *numOfVisited);
    for (int i = 0; i < *numOfVisited; ++i) {
      assert(temp[i] >= 0 && temp[i] < A->m);
      if (i > 0 && temp[i] == temp[i - 1]) {
        printf("%d duplicated\n", temp[i]);
        assert(false);
      }
    }

    delete[] temp;
  }
#endif

  return numLevels;
}

/**
 * Find minimum degree node among unvisited nodes.
 * Unvisited nodes are specified by color array
 */
static int getMinDegreeNode(const CSR *A, const int *nodes, int numOfNodes)
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
  for (int i = 0; i < numOfNodes; ++i) {
    int u = nodes[i];
    int degree = A->rowPtr[u + 1] - A->rowPtr[u];
    if (degree < temp_min) {
      temp_min = degree;
      temp_min_idx = u;
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

int getVerticesAtLevel(
  int *vertices,
  const int *levels, int level, const int *components, int numOfComponents)
{
  int idx = 0;
#pragma omp parallel for
  for (int i = 0; i < numOfComponents; ++i) {
    int u = components[i];
    if (levels[u] == level) {
      int idx_save = __sync_fetch_and_add(&idx, 1);
      vertices[idx_save] = u;
    }
  }

  sort(vertices, vertices + idx);

  return idx;
}

static void initializeLevels(
  int *levels, int m, const int *nodes, int numOfNodes)
{
#pragma omp parallel for
  for (int i = 0; i < numOfNodes; ++i) {
    assert(nodes[i] >= 0 && nodes[i] < m);
    levels[nodes[i]] = INT_MAX;
  }
}

int selectSourcesWithPseudoDiameter(
  const CSR *A, int *levels, const int *components, int numOfComponents, bfsAuxData *aux)
{
  // find the min degree node of this connected component
  int s = getMinDegreeNode(A, components, numOfComponents);
#ifdef PRINT_DBG
  printf("%d is the min degree node of this component\n", s);
#endif
  int e = -1;

  int *candidates = aux->candidates;

  while (e == -1) {
    initializeLevels(levels, A->m, components, numOfComponents);
    int diameter = bfs(A, s, levels, aux);

#ifdef PRINT_DBG
    printf("diameter from %d is %d\n", s, diameter);
#endif
    int nCandidates = getVerticesAtLevel(
      candidates, levels, diameter - 2, components, numOfComponents);

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

#ifdef PRINT_DBG
    printf("has %d candidates, brought down to %d\n", nCandidates, outIdx);
#endif

    nCandidates = outIdx;

    int minWidth = INT_MAX;
    for (int i = 0; i < nCandidates; ++i) {
      initializeLevels(levels, A->m, components, numOfComponents);

      int width = INT_MIN;
      int newDiameter = bfs(
        A, candidates[i], levels, aux, NULL, NULL, &width, &minWidth);
#ifdef PRINT_DBG
      printf("(diameter, width) from %d is (%d, %d)\n", candidates[i], newDiameter, width);
#endif
      if (-1 == newDiameter) { // short circuited
        continue;
      }
      else if (newDiameter > diameter && width < minWidth) {
        s = candidates[i];
#ifdef PRINT_DBG
        printf("select %d as the new starting point\n", s);
#endif
        e = -1;
        break;
      }
      else if (width < minWidth) {
        minWidth = width;
        e = candidates[i];
      }
    }
  } // iterate to find maximal diameter

  return s;
}

void CSR::getRCMPermutation(int *perm, int *inversePerm, int source /*=-1*/) const
{
  // 1. Start vertex
  double t;
  double sourceSelectionTime = 0, bfsTime = 0, prefixTime = 0, placeTime = 0;

  int *levels = new int[m];
  int maxDegree = INT_MIN;
#pragma omp parallel for reduction(max:maxDegree)
  for (int i = 0; i < m; ++i) {
    levels[i] = INT_MAX;
    maxDegree = max(maxDegree, rowPtr[i + 1] - rowPtr[i]);
  }
  int *children_array = new int[omp_get_max_threads()*maxDegree];

  bool sourcePreselected = source != -1;
  int *components = NULL;
  if (!sourcePreselected) components = new int[m];
  int offset = 0;
  int cnt = 0;

  bfsAuxData aux(m);
  volatile int *read_offset = new int[m + 1];
  volatile int *write_offset = new int[m + 1];
  int *prefixSum = new int[m + 1];

  // for each connected component
  int iEnd = sourcePreselected ? 1 : m;
  for (int i = 0; i < m; ++i) {
    // 1. Automatic selection of source
    int numOfComponents;
    if (!sourcePreselected) {
      if (levels[i] != INT_MAX) continue;
      ++cnt;

      t = omp_get_wtime();

      // collect nodes of this connected component
      bfs<true>(this, i, levels, &aux, components, &numOfComponents);
#ifdef PRINT_DBG
      printf("numOfComponents = %d\n", numOfComponents);
#endif

      // select source
      source = selectSourcesWithPseudoDiameter(
        this, levels, components, numOfComponents, &aux);
#ifdef PRINT_DBG
      printf("source selection takes %f\n", omp_get_wtime() - t);
      printf("source = %d\n", source);
#endif
      assert(source >= 0 && source < m);

      sourceSelectionTime += omp_get_wtime() - t;
    }
    else {
      numOfComponents = m;
    }

    // 2. BFS
    t = omp_get_wtime();
    if (!sourcePreselected) {
      initializeLevels(levels, m, components, numOfComponents);
    }
    int numLevels = bfs(this, source, levels, &aux);
#ifdef PRINT_DBG
    printf("numLevels = %d\n", numLevels);
    printf("bfs takes %f\n", omp_get_wtime() - t);
#endif
    bfsTime += omp_get_wtime() - t;

    // 3. Reorder
    t = omp_get_wtime();
    if (sourcePreselected) {
      prefixSumOfLevels(prefixSum, this, levels, numLevels, NULL, m);
    }
    else {
      prefixSumOfLevels(
        prefixSum, this, levels, numLevels, components, numOfComponents);
    }
#ifdef PRINT_DBG
    printf("prefix sum takes %f\n", omp_get_wtime() - t);
#endif
    prefixTime += omp_get_wtime() - t;

    t = omp_get_wtime();
    inversePerm[offset] = source;

    read_offset[0] = offset;
    write_offset[0] = offset + 1;

#pragma omp parallel
    {
      int nthreads = omp_get_num_threads();
      int tid = omp_get_thread_num();

      int *children = children_array + tid*maxDegree;

#pragma omp for
      for (int l = 1; l <= numLevels; ++l) {
        read_offset[l] = prefixSum[l] + offset;
        write_offset[l] = prefixSum[l] + offset;
      }

      for (int l = tid; l < numLevels; l += nthreads) {
        while (read_offset[l] != prefixSum[l + 1] + offset) {
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

#ifdef PRINT_DBG
    printf("place takes %f\n", omp_get_wtime() - t);
#endif
    placeTime += omp_get_wtime() - t;

    offset += numOfComponents;
  }

  delete[] levels;
  delete[] children_array;
  delete[] read_offset;
  delete[] write_offset;
  delete[] prefixSum;

  int *temp = new int[m];
#pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    temp[m - i - 1] = inversePerm[i];
  }
#pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    inversePerm[i] = temp[i];
    perm[inversePerm[i]] = i;
  }
  delete[] temp;

  printf("num of connected components = %d\n", cnt);
#ifdef PRINT_DBG
  printf("sourceSelectionTime = %f\n", sourceSelectionTime);
  printf("bfsTime = %f\n", bfsTime);
  printf("prefixTime = %f\n", prefixTime);
  printf("placeTime = %f\n", placeTime);
#endif
}
