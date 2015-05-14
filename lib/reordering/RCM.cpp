#include <cassert>
#include <climits>
#include <cstring>
#include <cstdio>

#include <vector>
#include <algorithm>

#include <omp.h>

#include "../CSR.hpp"
#include "../CSR_Interface.h"

using namespace std;

void getInversePerm(int *inversePerm, const int *perm, int n);

// compute prefix sum of levels
static void prefixSumOfLevels(
  int *prefixSum,
  const CSR *A, const int *levels, int numLevels,
  const int *components, int sizeOfComponents,
  bool parallel = true)
{
  if (parallel) {
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
        for (int i = 0; i < sizeOfComponents; ++i) {
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
        for (int t = 1; t < nthreads; ++t) {
          local_sum_array[t + 1] += local_sum_array[t];
        }
        assert(local_sum_array[nthreads] == sizeOfComponents);
        prefixSum[numLevels] = sizeOfComponents;
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
  } // parallel
  else {
    memset(prefixSum, 0, sizeof(int)*(numLevels + 1));

    if (NULL == components) {
      for (int i = 0; i < A->m; ++i) {
        assert(levels[i] != INT_MAX);
        ++prefixSum[levels[i] + 1];
      }
    }
    else {
      for (int i = 0; i < sizeOfComponents; ++i) {
        assert(components[i] >= 0 && components[i] < A->m);
        assert(levels[components[i]] != INT_MAX);
        ++prefixSum[levels[components[i]] + 1];
      }
    }

    for (int l = 0; l < numLevels; ++l) {
      prefixSum[l + 1] += prefixSum[l];
    }
  }
}

struct bfsAuxData
{
  int *q[2];
  int *qTail[2];
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

  qTail[0] = new int[omp_get_max_threads()];
  qTail[1] = new int[omp_get_max_threads()];

  qTailPrefixSum = new int[omp_get_max_threads() + 1];

  rowPtrs = new int[omp_get_max_threads()*m];
  nnzPrefixSum = new int[omp_get_max_threads() + 1];

  candidates = new int[omp_get_max_threads()*m];
}

bfsAuxData::~bfsAuxData()
{
  delete[] q[0];
  delete[] q[1];
  delete[] qTail[0];
  delete[] qTail[1];
  delete[] qTailPrefixSum;
  delete[] rowPtrs;
  delete[] nnzPrefixSum;
  delete[] candidates;
}

typedef char BitVectorType;

class BitVector
{
public :
  BitVector(int n) {
    int n_ = (n + sizeof(BitVectorType) - 1)/sizeof(BitVectorType);
    bv_ = new BitVectorType[n_];
#pragma omp parallel for
    for (int i = 0; i < n_; ++i) {
      bv_[i] = 0;
    }
  }

  ~BitVector() { delete[] bv_; }

  void set(int i) {
    bv_[getIndexOf_(i)] |= getMaskOf_(i);
  }

  void atomicSet(int i) {
    BitVectorType mask = getMaskOf_(i);
    __sync_fetch_and_or(bv_ + getIndexOf_(i), mask);
  }

  bool get(int i) const {
    return bv_[getIndexOf_(i)] & getMaskOf_(i);
  }

  bool testAndSet(int i) {
    if (!get(i)) {
      BitVectorType mask = getMaskOf_(i);
      BitVectorType prev = __sync_fetch_and_or(bv_ + getIndexOf_(i), mask);
      return !(prev & mask);
    }
    else {
      return false;
    }
  }

  bool atomicClear(int i) {
    __sync_fetch_and_and(bv_ + getIndexOf_(i), ~getMaskOf_(i));
  }

private :
  typedef char BitVectorType;

  static int getIndexOf_(int i) { return i/sizeof(BitVectorType); }
  static int getBitIndexOf_(int i) { return i%sizeof(BitVectorType); }
  static BitVectorType getMaskOf_(int i) { return 1 << getBitIndexOf_(i); }

  BitVectorType *bv_;
};

/**
 * @return -1 if shortcircuited num of levels otherwise
 *
 * pre-condition: levels should be initialized to -1
 */
template<bool SET_LEVEL = true, bool OUTPUT_VISITED = false>
int bfs_serial(
  const CSR *A, int source, int *levels, BitVector *bv,
  bfsAuxData *aux,
  int *visited = NULL, int *numOfVisited = NULL,
  int *width = NULL, int *shortCircuitWidth = NULL) {

  int tid = omp_get_thread_num();

  int numLevels = 0;
  if (SET_LEVEL) levels[source] = numLevels;
  bv->atomicSet(source);

  int **q = aux->q;
  q[0][tid*A->m] = source;

  int *qTail[2] = { aux->qTail[0], aux->qTail[1] };
  qTail[0][tid] = 1;

  if (OUTPUT_VISITED) *numOfVisited = 0;

  while (true) {
    if (width) {
      *width = max(*width, qTail[numLevels%2][tid]);
    }
    if (shortCircuitWidth) {
      if (qTail[numLevels%2][tid] > *shortCircuitWidth) {
        numLevels = -1;
      }
    }
    if (OUTPUT_VISITED) *numOfVisited += qTail[numLevels%2][tid];

    ++numLevels;

    if (qTail[1 - numLevels%2][tid] == 0 || numLevels == -1) break;

    if (OUTPUT_VISITED) {
      memcpy(
        visited + *numOfVisited - qTail[1 - numLevels%2][tid],
        q[1 - numLevels%2] + tid*A->m,
        sizeof(int)*qTail[1 - numLevels%2][tid]);
    }

    int *tailPtr = q[numLevels%2] + tid*A->m;

    for (int i = 0; i < qTail[1 - numLevels%2][tid]; ++i) {
      int u = q[1 - numLevels%2][i + tid*A->m];
      assert(!SET_LEVEL || levels[u] == numLevels - 1);

      for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
        int v = A->colIdx[j];
        if (!bv->get(v)) {
          bv->atomicSet(v);
          if (SET_LEVEL) levels[v] = numLevels;

          *tailPtr = v;
          ++tailPtr;
        }
      }
    } // for each current node u

    qTail[numLevels%2][tid] = tailPtr - (q[numLevels%2] + tid*A->m);
  } // while true

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
 * @return -1 if shortcircuited num of levels otherwise
 *
 * pre-condition: levels should be initialized to -1
 */
template<bool SET_LEVEL = true, bool OUTPUT_VISITED = false>
int bfs(
  const CSR *A, int source, int *levels, BitVector *bv,
  bfsAuxData *aux,
  int *visited = NULL, int *numOfVisited = NULL,
  int *width = NULL, int *shortCircuitWidth = NULL) {

  int numLevels = 0;
  if (SET_LEVEL) levels[source] = numLevels;
  bv->set(source);

  int **q = aux->q;
  q[0][0] = source;

  int *qTail[2] = { aux->qTail[0], aux->qTail[1] };
  qTail[0][0] = 1;

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
    qTail[0][tid] = 0;
    nnzPrefixSum[tid + 1] = 0;
  }

  while (true) {
#pragma omp barrier
#pragma omp master
    {
      for (int t = 0; t < nthreads; ++t) {
        qTailPrefixSum[t + 1] = qTailPrefixSum[t] + qTail[numLevels%2][t];
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
      iBegin = qTail[1 - numLevels%2][tEnd - 1];
    }
    else {
      iBegin = upper_bound(
          rowPtrs + tBegin*A->m, rowPtrs + tBegin*A->m + qTail[1 - numLevels%2][tBegin],
          nnzPerThread*tid - nnzPrefixSum[tBegin]) -
        (rowPtrs + tBegin*A->m) - 1;
    }

    if (tEnd == nthreads) {
      iEnd = 0;
    }
    else {
      iEnd = upper_bound(
          rowPtrs + tEnd*A->m, rowPtrs + tEnd*A->m + qTail[1 - numLevels%2][tEnd],
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
          i < (t == tEnd ? iEnd : qTail[1 - numLevels%2][t]);
          ++i) {
        int u = q[1 - numLevels%2][t*A->m + i];
        assert(!SET_LEVEL || levels[u] == numLevels - 1);

        for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
          int v = A->colIdx[j];
          if (OUTPUT_VISITED) {
            if (bv->testAndSet(v)) {
              if (SET_LEVEL) levels[v] = numLevels;

              *tailPtr = v;
              *(rowPtr + 1) = *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

              ++tailPtr;
              ++rowPtr;
            }
          }
          else {
            if (!bv->get(v)) {
              bv->set(v);
              if (SET_LEVEL) levels[v] = numLevels;

              *tailPtr = v;
              *(rowPtr + 1) = *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

              ++tailPtr;
              ++rowPtr;
            }
          }
        }
      } // for each current node u
    }

    qTail[numLevels%2][tid] = tailPtr - (q[numLevels%2] + tid*A->m);
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
static int getMinDegreeNode(const CSR *A, const int *nodes, int numOfNodes, bool parallel = true)
{
  int global_min_idx;

  if (parallel) {
    int local_min[omp_get_max_threads()];
    int local_min_idx[omp_get_max_threads()];

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
  }
  else {
    int global_min = INT_MAX;
    global_min_idx = -1;

    for (int i = 0; i < numOfNodes; ++i) {
      int u = nodes[i];
      int degree = A->rowPtr[u + 1] - A->rowPtr[u];
      if (degree < global_min) {
        global_min = degree;
        global_min_idx = u;
      }
    }
  }

  return global_min_idx;
}

static void initializeBitVector(
  BitVector *bv, int m, const int *nodes, int numOfNodes, bool parallel = true)
{
  if (parallel) {
#pragma omp parallel for
    for (int i = 0; i < numOfNodes; ++i) {
      int u = nodes[i];
      assert(u >= 0 && u < m);
      bv->atomicClear(u);
    }
  }
  else {
    for (int i = 0; i < numOfNodes; ++i) {
      int u = nodes[i];
      assert(u >= 0 && u < m);
      bv->atomicClear(u);
    }
  }
}

int selectSourcesWithPseudoDiameter(
  const CSR *A, BitVector *bv, const int *components, int sizeOfComponents, bfsAuxData *aux)
{
  // find the min degree node of this connected component
  int s = getMinDegreeNode(A, components, sizeOfComponents);
//#define PRINT_DBG
#ifdef PRINT_DBG
  printf("%d is the min degree node of this component\n", s);
#endif
  int e = -1;

  int tid = omp_get_thread_num();
  int *candidates = aux->candidates;

  while (e == -1) {
    initializeBitVector(bv, A->m, components, sizeOfComponents);
    int diameter = bfs<false>(A, s, NULL, bv, aux);

#ifdef PRINT_DBG
    printf("diameter from %d is %d\n", s, diameter);
#endif

    int nCandidates = 0;
    for (int t = 0; t < omp_get_max_threads(); ++t) {
      for (int j = 0; j < aux->qTail[diameter%2][t]; ++j) {
        int u = aux->q[diameter%2][t*A->m + j];
        candidates[nCandidates] = u;
        ++nCandidates;
      }
    }

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
      for (int k = 0; !adjacent && k < outIdx; ++k) {
        if (candidates[k] == u) adjacent = true;
      }
      for (int j = A->rowPtr[u]; !adjacent && j < A->rowPtr[u + 1]; ++j) {
        int v = A->colIdx[j];
        for (int k = 0; !adjacent && k < outIdx; ++k) {
          if (candidates[k] == v) adjacent = true;
        }
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
      initializeBitVector(bv, A->m, components, sizeOfComponents);

      int width = INT_MIN;
      int newDiameter =
        bfs<false>(
          A, candidates[i], NULL, bv, aux, NULL, NULL, &width, &minWidth);

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

class DegreeComparator
{
public :
  DegreeComparator(const int *rowPtr) : rowPtr_(rowPtr) { };

  bool operator()(int a, int b) {
    return rowPtr_[a + 1] - rowPtr_[a] < rowPtr_[b + 1] - rowPtr_[b];
  }

private :
  const int *rowPtr_;
};

void CSR::getRCMPermutation(int *perm, int *inversePerm, int source /*=-1*/)
{
  assert(isSymmetric(false)); // check structural symmetry

  int oldBase = base;
  make0BasedIndexing();

  // 1. Start vertex
  double t;
  double sourceSelectionTime1 = 0, bfsTime1 = 0, prefixTime1 = 0, placeTime1 = 0;
  double sourceSelectionTime2 = 0, bfsTime2 = 0, prefixTime2 = 0, placeTime2 = 0;

  int *levels = new int[m];
  int maxDegree = INT_MIN;
#pragma omp parallel for reduction(max:maxDegree)
  for (int i = 0; i < m; ++i) {
    maxDegree = max(maxDegree, rowPtr[i + 1] - rowPtr[i]);
  }
#ifdef PRINT_DBG
  printf("maxDegree = %d\n", maxDegree);
#endif

  BitVector *bv = new BitVector(m);

  int *children_array = new int[omp_get_max_threads()*maxDegree];

  int *components_array = new int[omp_get_max_threads()*m];
  int singletonCnt = 0, twinCnt = 0;

  bfsAuxData aux(m);
  volatile int *write_offset = new int[(m + 1)*16];
  int *prefixSum_array = new int[omp_get_max_threads()*(m + 1)];

  DegreeComparator comparator(rowPtr);

  int numOfComponents;
  int *compToRoot, *compSizes, *compSizePrefixSum;

  CSR_FindConnectedComponents(
    (const CSR_Handle *)this,
    &numOfComponents, &compToRoot, &compSizes, &compSizePrefixSum);

  const int PAR_THR = 16;

#define MEASURE_LOAD_BALANCE
#ifdef MEASURE_LOAD_BALANCE
  double barrierTimes[omp_get_max_threads()];
  double tBegin = omp_get_wtime();
  double barrierTimeSum = 0;
#endif

#pragma omp parallel reduction(+:singletonCnt,twinCnt,sourceSelectionTime1,bfsTime1,prefixTime1,placeTime1)
  {
    int tid = omp_get_thread_num();

    int *children = children_array + tid*maxDegree;
    int *components = components_array + tid*m;
    int *prefixSum = prefixSum_array + tid*(m + 1);

  // for each connected component
#pragma omp for
  for (int c = 0; c < numOfComponents; ++c) {
    int i = compToRoot[c];
    int offset = compSizePrefixSum[c];

    // 1. Automatic selection of source
    int sizeOfComponents;

    if (compSizes[c] >= PAR_THR) continue;

    // short circuit for a singleton or a twin
    if (rowPtr[i + 1] == rowPtr[i] + 1 && colIdx[rowPtr[i]] == i || rowPtr[i + 1] == rowPtr[i]) {
      inversePerm[m - offset - 1] = i;
      perm[i] = m - offset - 1;
      ++singletonCnt;
      continue;
    }
    else if (rowPtr[i + 1] == rowPtr[i] + 1) {
      int u = colIdx[rowPtr[i]];
      if (rowPtr[u + 1] == rowPtr[u] + 1 && colIdx[rowPtr[u]] == i) {
        inversePerm[m - offset - 1] = i;
        inversePerm[m - (offset + 1) - 1] = u;
        perm[i] = m - offset - 1;
        perm[u] = m - (offset + 1) - 1;
        bv->atomicSet(u);
        ++twinCnt;
        continue;
      }
    }
    else if (rowPtr[i + 1] == rowPtr[i] + 2) {
      int u = -1;
      if (colIdx[rowPtr[i]] == i) {
        u = colIdx[rowPtr[i] + 1];
      }
      else if (colIdx[rowPtr[i] + 1] == i) {
        u = colIdx[rowPtr[i]];
      }
      if (u != -1 &&
        rowPtr[u + 1] == rowPtr[u] + 2 &&
          (colIdx[rowPtr[u]] == u && colIdx[rowPtr[u] + 1] == i ||
            colIdx[rowPtr[u] + 1] == u && colIdx[rowPtr[u]] == i)) {
        inversePerm[m - offset - 1] = i;
        inversePerm[m - (offset + 1) - 1] = u;
        perm[i] = m - offset - 1;
        perm[u] = m - (offset + 1) - 1;
        bv->atomicSet(u);
        ++twinCnt;
        continue;
      }
    }

    double t = omp_get_wtime();

    // collect nodes of this connected component
    bfs_serial<false, true>(this, i, NULL, bv, &aux, components, &sizeOfComponents);
    assert(sizeOfComponents == compSizes[c]);
    if (sizeOfComponents != compSizes[c])
    {
      printf("sizeOfComponent %d expected %d actual\n", compSizes[c], sizeOfComponents);
    }
#ifdef PRINT_DBG
    printf("sizeOfComponents = %d\n", sizeOfComponents);
#endif

    // select source
    int source = i;
    sourceSelectionTime1 += omp_get_wtime() - t;

    // 2. BFS
    t = omp_get_wtime();
    initializeBitVector(bv, m, components, sizeOfComponents, false);
    int numLevels = bfs_serial(this, source, levels, bv, &aux);
#ifdef PRINT_DBG
    printf("numLevels = %d\n", numLevels);
    printf("bfs takes %f\n", omp_get_wtime() - t);
#endif
    bfsTime1 += omp_get_wtime() - t;

    // 3. Reorder
    t = omp_get_wtime();
    prefixSumOfLevels(
      prefixSum, this, levels, numLevels, components, sizeOfComponents, false);
#ifdef PRINT_DBG
    printf("prefix sum takes %f\n", omp_get_wtime() - t);
#endif
    prefixTime1 += omp_get_wtime() - t;

    t = omp_get_wtime();
    inversePerm[m - offset - 1] = source;
    perm[source] = m - offset - 1;

    for (int l = 0; l < numLevels; ++l) {
      int r = prefixSum[l] + offset;
      int w = prefixSum[l + 1] + offset;
      while (r != prefixSum[l + 1] + offset) {
        int u = inversePerm[m - r - 1];
        ++r;
        int childrenIdx = 0;
        for (int j = rowPtr[u]; j < rowPtr[u + 1]; ++j) {
          int v = colIdx[j];
          if (levels[v] == l + 1) {
            children[childrenIdx] = v;
            ++childrenIdx;
            levels[v] = -1;
          }
        }

        std::sort(children, children + childrenIdx, comparator);

        for (int i = 0; i < childrenIdx; ++i) {
          int c = children[i];
          int idx = m - (w + i) - 1;
          inversePerm[idx] = c;
          perm[c] = idx;
        }
        w += childrenIdx;
      }
    }

#ifdef PRINT_DBG
    printf("place takes %f\n", omp_get_wtime() - t);
#endif
    placeTime1 += omp_get_wtime() - t;
  }

#ifdef MEASURE_LOAD_BALANCE
    double t = omp_get_wtime();
#pragma omp barrier
    barrierTimes[tid] = omp_get_wtime() - t;

#pragma omp master
    {
      double tEnd = omp_get_wtime();
      for (int i = 0; i < omp_get_num_threads(); ++i) {
        barrierTimeSum += barrierTimes[i];
      }
      printf("%f load imbalance = %f\n", tEnd - tBegin, barrierTimeSum/(tEnd - tBegin)/omp_get_num_threads());
    }
#undef MEASURE_LOAD_BALANCE
#endif // MEASURE_LOAD_BALANCE
  } // omp parallel

  double timeFirstPhase = omp_get_wtime() - tBegin;
  tBegin = omp_get_wtime();

  int *components = components_array;
  int *prefixSum = prefixSum_array;

  // for each connected component
  for (int c = 0; c < numOfComponents; ++c) {
    int i = compToRoot[c];
    int offset = compSizePrefixSum[c];

    // 1. Automatic selection of source
    int sizeOfComponents;

    if (compSizes[c] < PAR_THR) continue;

    t = omp_get_wtime();

    // collect nodes of this connected component
    bfs<false, true>(this, i, NULL, bv, &aux, components, &sizeOfComponents);
    assert(sizeOfComponents == compSizes[c]);
    if (sizeOfComponents != compSizes[c])
    {
      printf("sizeOfComponent %d expected %d actual\n", compSizes[c], sizeOfComponents);
    }
#ifdef PRINT_DBG
    printf("sizeOfComponents = %d\n", sizeOfComponents);
#endif

    // select source
    source = selectSourcesWithPseudoDiameter(
      this, bv, components, sizeOfComponents, &aux);
#ifdef PRINT_DBG
    printf("source selection takes %f\n", omp_get_wtime() - t);
    printf("source = %d\n", source);
#endif
    assert(source >= 0 && source < m);

    sourceSelectionTime2 += omp_get_wtime() - t;

    // 2. BFS
    t = omp_get_wtime();
    initializeBitVector(bv, m, components, sizeOfComponents);
    int numLevels = bfs(this, source, levels, bv, &aux);
#ifdef PRINT_DBG
    printf("numLevels = %d\n", numLevels);
    printf("bfs takes %f\n", omp_get_wtime() - t);
#endif
    bfsTime2 += omp_get_wtime() - t;

    // 3. Reorder
    t = omp_get_wtime();
    prefixSumOfLevels(
      prefixSum, this, levels, numLevels, components, sizeOfComponents);
#ifdef PRINT_DBG
    printf("prefix sum takes %f\n", omp_get_wtime() - t);
#endif
    prefixTime2 += omp_get_wtime() - t;

    t = omp_get_wtime();
    inversePerm[m - offset - 1] = source;
    perm[source] = m - offset - 1;

#pragma omp parallel
    {
      int nthreads = omp_get_num_threads();
      int tid = omp_get_thread_num();

      int *children = children_array + tid*maxDegree;

      for (int l = tid; l <= numLevels; l += nthreads) {
        write_offset[16*l] = prefixSum[l] + offset;
      }
      if (0 == tid) {
        write_offset[0] = offset + 1;
      }

#pragma omp barrier

      for (int l = tid; l < numLevels; l += nthreads) {
        int r = prefixSum[l] + offset;
        while (r != prefixSum[l + 1] + offset) {
          while (r == write_offset[16*l]); // spin
          int u = inversePerm[m - r - 1];
          ++r;
          int childrenIdx = 0;
          for (int j = rowPtr[u]; j < rowPtr[u + 1]; ++j) {
            int v = colIdx[j];
            if (levels[v] == l + 1) {
              children[childrenIdx] = v;
              ++childrenIdx;
              levels[v] = -1;
            }
          }

          std::sort(children, children + childrenIdx, comparator);

          int w = write_offset[16*(l + 1)];
          for (int i = 0; i < childrenIdx; ++i) {
            int c = children[i];
            int idx = m - (w + i) - 1;
            inversePerm[idx] = c;
            perm[c] = idx;
            write_offset[16*(l + 1)] = w + i + 1;
          }
        }
      } // for each level
    } // omp parallel

#ifdef PRINT_DBG
    printf("place takes %f\n", omp_get_wtime() - t);
#endif
    placeTime2 += omp_get_wtime() - t;
  }

  double timeSecondPhase = omp_get_wtime() - tBegin;

  delete[] levels;
  delete bv;
  delete[] children_array;
  delete[] components;
  delete[] write_offset;
  delete[] prefixSum;

  printf("num of connected components = %d (singleton = %d, twin = %d)\n", numOfComponents, singletonCnt, twinCnt);
  printf("firstPhaseTime (parallel over components) = %f\n", timeFirstPhase);
  printf("\tsourceSelectionTime = %f\n", sourceSelectionTime1/omp_get_max_threads());
  printf("\tbfsTime = %f\n", bfsTime1/omp_get_max_threads());
  printf("\tprefixTime = %f\n", prefixTime1/omp_get_max_threads());
  printf("\tplaceTime = %f\n", placeTime1/omp_get_max_threads());
  printf("\tloadImbalanceTime = %f\n", barrierTimeSum/omp_get_max_threads());
  printf("secondPhaseTime (parallel within components) = %f\n", timeSecondPhase);
  printf("\tsourceSelectionTime = %f\n", sourceSelectionTime2);
  printf("\tbfsTime = %f\n", bfsTime2);
  printf("\tprefixTime = %f\n", prefixTime2);
  printf("\tplaceTime = %f\n", placeTime2);

  if (1 == oldBase) {
    make1BasedIndexing();
  }
}
