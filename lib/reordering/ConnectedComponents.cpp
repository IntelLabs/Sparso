#include <cassert>
#include <algorithm>

#include <omp.h>

#include <CSR.hpp>
#include <CSR_Interface.h>

using namespace std;

template<int BASE = 0>
void findConnectedComponents_(
  const CSR *A,
  int *numOfComponents, int **compToRoot, int **compSizes, int **compSizePrefixSum)
{
  volatile int *p = new int[A->m];

  int cnts[omp_get_max_threads() + 1]; // prefix sum of # of connected components
  cnts[0] = 0;

  *compToRoot = NULL;
  int *rootToComp = new int[A->m];
  *compSizes = NULL;
  *numOfComponents = 0;

  double t = omp_get_wtime();

#pragma omp parallel
  {
  int tid = omp_get_thread_num();
  int nthreads = omp_get_num_threads();

  int iPerThread = (A->m + nthreads - 1)/nthreads;
  int iBegin = min(iPerThread*tid, A->m);
  int iEnd = min(iBegin + iPerThread, A->m);

  for (int i = iBegin; i< iEnd; ++i) {
    p[i] = i;
  }

#pragma omp barrier

  int nnz = A->rowPtr[A->m] - BASE;
  int nnzPerThread = (nnz + nthreads - 1)/nthreads;
  int xBegin = lower_bound(A->rowPtr, A->rowPtr + A->m, nnzPerThread*tid + BASE) - A->rowPtr;
  int xEnd = lower_bound(A->rowPtr, A->rowPtr + A->m, nnzPerThread*(tid + 1) + BASE) - A->rowPtr;
  assert(xBegin <= xEnd);
  assert(xBegin >= 0 && xBegin <= A->m);
  assert(xEnd >= 0 && xEnd <= A->m);

  for (int x = xBegin; x < xEnd; ++x) {
    for (int j = A->rowPtr[x] - BASE; j < A->rowPtr[x + 1] - BASE; ++j) {
      int y = A->colIdx[j];
      if (p[x] != p[y]) {
        // union
        int r_x = x, r_y = y;
        while (true) {
          int old_p_r_x = p[r_x]; int old_p_r_y = p[r_y];
          if (old_p_r_x == old_p_r_y) break;

          int old_r_x = r_x; int old_r_y = r_y;

          r_x = old_p_r_x > old_p_r_y ? old_r_x : old_r_y;
          r_y = old_p_r_x > old_p_r_y ? old_r_y : old_r_x;
          int p_r_x = old_p_r_x > old_p_r_y ? old_p_r_x : old_p_r_y;
          int p_r_y = old_p_r_x > old_p_r_y ? old_p_r_y : old_p_r_x;

          if (p_r_x == r_x && __sync_bool_compare_and_swap(&p[r_x], r_x, p_r_y)) {
            break;
          }
          p[r_x] = p_r_y;
          r_x = p_r_x;
        } // while
      } // p[x] != p[y]
    }
  } // for each row x

#pragma omp barrier

  // path compression so that all p[i] points to its root
  // and count # of components
  int compId = 0;
  for (int i = iBegin; i < iEnd; ++i) {
    int r = i;
    while (p[r] != r) {
      r = p[r];
    }
    p[i] = r;
    if (r == i) ++compId;
  }

  cnts[tid + 1] = compId;

  // prefix sum # of components
#pragma omp barrier
#pragma omp master
  {
    for (int i = 1; i < nthreads; ++i) {
      cnts[i + 1] += cnts[i];
    }
    *numOfComponents = cnts[nthreads];
    *compToRoot = new int[*numOfComponents];
    *compSizes = new int[*numOfComponents*nthreads];
    *compSizePrefixSum = new int[*numOfComponents];
  }
#pragma omp barrier

  // compId <-> root map
  compId = cnts[tid];
  for (int i = iBegin; i < iEnd; ++i) {
    int r = p[i];
    if (r == i) {
      (*compToRoot)[compId] = r;
      rootToComp[r] = compId;
      ++compId;
    }
  }
  
#pragma omp barrier

  // count thread-private component sizes
  for (int c = 0; c < *numOfComponents; ++c) {
    (*compSizes)[c + tid*(*numOfComponents)] = 0;
  }

  for (int i = iBegin; i < iEnd; ++i) {
    (*compSizes)[rootToComp[p[i]] + tid*(*numOfComponents)]++;
  }

#pragma omp barrier

  int cPerThread = (*numOfComponents + nthreads - 1)/nthreads;
  int cBegin = min(cPerThread*tid, *numOfComponents);
  int cEnd = min(cBegin + cPerThread, *numOfComponents);

  // reduce component sizes
  int localCnt = 0;
  for (int c = cBegin; c < cEnd; ++c) {
    int sum = (*compSizes)[c];
    for (int t = 1; t < nthreads; ++t) {
      sum += (*compSizes)[c + t*(*numOfComponents)];
    }
    (*compSizes)[c] = sum;
    localCnt += sum;
  }
  cnts[tid + 1] = localCnt;

#pragma omp barrier
#pragma omp master
  {
    for (int i = 1; i < nthreads - 1; ++i) {
      cnts[i + 1] += cnts[i];
    }
  }
#pragma omp barrier

  localCnt = cnts[tid];
  for (int c = cBegin; c < cEnd; ++c) {
    (*compSizePrefixSum)[c] = localCnt;
    localCnt += (*compSizes)[c];
  }
  } // omp parallel

  printf("finding connected components takes %f\n", omp_get_wtime() - t);

#ifndef NDEBUG
  int cnt = 0;
  for (int i = 0; i < (*numOfComponents); ++i) {
    cnt += (*compSizes)[i];
  }
  assert(cnt == A->m);
#endif

  printf("num of connected components = %d\n", (*numOfComponents));
}

extern "C" {

void CSR_FindConnectedComponents(
  const CSR_Handle *A,
  int *numOfComponents, int **compToRoot, int **compSizes, int **compSizePrefixSum)
{
  if (0 == ((CSR *)A)->base) {
    findConnectedComponents_<0>((const CSR *)A, numOfComponents, compToRoot, compSizes, compSizePrefixSum);
  }
  else {
    assert(1 == ((CSR *)A)->base);
    findConnectedComponents_<1>((const CSR *)A, numOfComponents, compToRoot, compSizes, compSizePrefixSum);
  }
}

} // extern "C"
