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
  const int *components, int sizeOfComponents)
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

  candidates = new int[m];
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
template<bool OUTPUT_VISITED = false>
int bfs(
  const CSR *A, int source, int *levels, BitVector *bv,
  bfsAuxData *aux,
  int *visited = NULL, int *numOfVisited = NULL,
  int *width = NULL, int *shortCircuitWidth = NULL) {

  int numLevels = 0;
  levels[source] = numLevels;
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
      iBegin = qTail[numLevels%2][tEnd - 1];
    }
    else {
      iBegin = upper_bound(
          rowPtrs + tBegin*A->m, rowPtrs + tBegin*A->m + qTail[numLevels%2][tBegin],
          nnzPerThread*tid - nnzPrefixSum[tBegin]) -
        (rowPtrs + tBegin*A->m) - 1;
    }

    if (tEnd == nthreads) {
      iEnd = 0;
    }
    else {
      iEnd = upper_bound(
          rowPtrs + tEnd*A->m, rowPtrs + tEnd*A->m + qTail[numLevels%2][tEnd],
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
          i < (t == tEnd ? iEnd : qTail[numLevels%2][t]);
          ++i) {
        int u = q[1 - numLevels%2][t*A->m + i];
        assert(levels[u] == numLevels - 1);

        for (int j = A->rowPtr[u]; j < A->rowPtr[u + 1]; ++j) {
          int v = A->colIdx[j];
          if (OUTPUT_VISITED) {
            if (bv->testAndSet(v)) {
              levels[v] = numLevels;

              *tailPtr = v;
              *(rowPtr + 1) = *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

              ++tailPtr;
              ++rowPtr;
            }
          }
          else {
            if (!bv->get(v)) {
              bv->set(v);
              if (INT_MAX == levels[v]) {
                levels[v] = numLevels;

                *tailPtr = v;
                *(rowPtr + 1) = *rowPtr + A->rowPtr[v + 1] - A->rowPtr[v];

                ++tailPtr;
                ++rowPtr;
              }
            }
          }
        }
      } // for each current node u
    }

#pragma omp barrier

    qTail[1 - numLevels%2][tid] = tailPtr - (q[numLevels%2] + tid*A->m);
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
  const int *levels, int level, const int *components, int sizeOfComponents)
{
  int idx = 0;
#pragma omp parallel for
  for (int i = 0; i < sizeOfComponents; ++i) {
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
  int *levels, BitVector *bv, int m, const int *nodes, int numOfNodes)
{
#pragma omp parallel for
  for (int i = 0; i < numOfNodes; ++i) {
    int u = nodes[i];
    assert(u >= 0 && u < m);
    levels[u] = INT_MAX;
    bv->atomicClear(u);
  }
}

int selectSourcesWithPseudoDiameter(
  const CSR *A, int *levels, BitVector *bv, const int *components, int sizeOfComponents, bfsAuxData *aux)
{
  // find the min degree node of this connected component
  int s = getMinDegreeNode(A, components, sizeOfComponents);
#define PRINT_DBG
#ifdef PRINT_DBG
  printf("%d is the min degree node of this component\n", s);
#endif
  int e = -1;

  int *candidates = aux->candidates;

  while (e == -1) {
    initializeLevels(levels, bv, A->m, components, sizeOfComponents);
    int diameter = bfs(A, s, levels, bv, aux);

#ifdef PRINT_DBG
    printf("diameter from %d is %d\n", s, diameter);
#endif
    int nCandidates = getVerticesAtLevel(
      candidates, levels, diameter - 2, components, sizeOfComponents);

    /*for (int i = 0; i < omp_get_max_threads(); ++i) {
      for (int j = 0; j < A->m; ++j) {
        int u = aux->q[diameter%2][i*A->m + j];
        if (u < 0 || u >= A->m || levels[u] != diameter - 2) {
          break;
        }
        printf("%d ", u);
      }
    }
    printf("\n");*/

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
      initializeLevels(levels, bv, A->m, components, sizeOfComponents);

      int width = INT_MIN;
      int newDiameter = bfs(
        A, candidates[i], levels, bv, aux, NULL, NULL, &width, &minWidth);
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
#ifdef PRINT_DBG
  printf("maxDegree = %d\n", maxDegree);
#endif

  BitVector *bv = new BitVector(m);

  int *children_array = new int[omp_get_max_threads()*maxDegree];

  bool sourcePreselected = source != -1;
  int *components = NULL;
  if (!sourcePreselected) components = new int[m];
  int offset = 0;
  int connectedCmpCnt = 0, singletonCnt = 0, twinCnt = 0;

  bfsAuxData aux(m);
  volatile int *write_offset = new int[(m + 1)*16];
  int *prefixSum = new int[m + 1];

  DegreeComparator comparator(rowPtr);

  // for each connected component
  int iEnd = sourcePreselected ? 1 : m;
  for (int i = 0; i < m; ++i) {
    // 1. Automatic selection of source
    int sizeOfComponents;
    if (!sourcePreselected) {

      if (bv->get(i)) continue;
      ++connectedCmpCnt;

      // short circuit for a singleton or a twin
      if (rowPtr[i + 1] == rowPtr[i] + 1 && colIdx[rowPtr[i]] == i || rowPtr[i + 1] == rowPtr[i]) {
        inversePerm[m - offset - 1] = i;
        perm[i] = m - offset - 1;
        ++offset;
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
          levels[u] = 1;
          bv->set(u);
          offset += 2;
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
          levels[u] = 1;
          bv->set(u);
          offset += 2;
          ++twinCnt;
          continue;
        }
      }

      t = omp_get_wtime();

      // collect nodes of this connected component
      bfs<true>(this, i, levels, bv, &aux, components, &sizeOfComponents);
#ifdef PRINT_DBG
      printf("sizeOfComponents = %d\n", sizeOfComponents);
#endif

      // select source
      source = selectSourcesWithPseudoDiameter(
        this, levels, bv, components, sizeOfComponents, &aux);
#ifdef PRINT_DBG
      printf("source selection takes %f\n", omp_get_wtime() - t);
      printf("source = %d\n", source);
#endif
      assert(source >= 0 && source < m);

      sourceSelectionTime += omp_get_wtime() - t;
    }
    else {
      sizeOfComponents = m;
    }

    // 2. BFS
    t = omp_get_wtime();
    if (!sourcePreselected) {
      initializeLevels(levels, bv, m, components, sizeOfComponents);
    }
    int numLevels = bfs(this, source, levels, bv, &aux);
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
        prefixSum, this, levels, numLevels, components, sizeOfComponents);
    }
#ifdef PRINT_DBG
    printf("prefix sum takes %f\n", omp_get_wtime() - t);
#endif
    prefixTime += omp_get_wtime() - t;

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
    placeTime += omp_get_wtime() - t;

    offset += sizeOfComponents;
  }

  delete[] levels;
  delete bv;
  delete[] children_array;
  delete[] write_offset;
  delete[] prefixSum;

  printf("num of connected components = %d (singleton = %d, twin = %d)\n", connectedCmpCnt, singletonCnt, twinCnt);
  printf("sourceSelectionTime = %f\n", sourceSelectionTime);
  printf("bfsTime = %f\n", bfsTime);
  printf("prefixTime = %f\n", prefixTime);
  printf("placeTime = %f\n", placeTime);
}
