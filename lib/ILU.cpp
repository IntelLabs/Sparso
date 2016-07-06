/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <cstring>
#include "ILU.hpp"
#include "SpMP/synk/barrier.hpp"

using namespace std;

namespace SpMP
{

void ilu0(double *lu, const CSR& A, const LevelSchedule& schedule)
{
  int base = A.getBase();

  const int *rowptr = A.rowptr - base;
  const int *colidx = A.colidx - base;
  const int *diagptr = A.diagptr - base;
  const double *values = A.values - base;

  lu -= base;

#pragma omp parallel
  {
    int tid = omp_get_thread_num();

#pragma omp for
    for (int i = base; i < A.getNnz() + base; i++) {
      lu[i] = values[i];
    }

    const int ntasks = schedule.ntasks;
    const short *nparents = schedule.nparentsForward;
    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    const int *perm = schedule.threadContToOrigPerm;

    int nBegin, nEnd;
    getSimpleThreadPartition(&nBegin, &nEnd, ntasks);

    volatile int *taskFinished = schedule.taskFinished;
    int **parents = schedule.parentsForward;

    memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

    synk::Barrier::getInstance()->wait(tid);

    for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
      SPMP_LEVEL_SCHEDULE_WAIT;

      for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
        int row = perm[i] + base;

        for (int j = rowptr[row]; j < diagptr[row]; ++j) {
          int c = colidx[j];
          double tmp = lu[j] /= lu[diagptr[c]];

          int k1 = j + 1, k2 = diagptr[c] + 1;

          while (k1 < rowptr[row + 1] && k2 < rowptr[c + 1]) {
            if (colidx[k1] < colidx[k2]) ++k1;
            else if (colidx[k1] > colidx[k2]) ++k2;
            else {
              lu[k1] -= tmp*lu[k2];
              ++k1; ++k2;
            }
          }
        }
      } // for each row

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each level
  } // omp parallel
}

/*void ilu0(CSR& L, CSR& U, const CSR& A, const LevelSchedule& schedule)
{
  int base = A.getBase();

  const int *rowptr = A.rowptr - base;
  const int *colidx = A.colidx - base;
  const int *diagptr = A.diagptr - base;
  const double *values = A.values - base;

  lu -= base;

  L.dealloc();
  U.dealloc();

  L.m = U.m = A.m;
  L.n = U.m = A.n;

  L.rowptr = MALLOC(int, L.m + 1);
  U.rowptr = MALLOC(int, L.m + 1);
  L.idiag = MALLOC(double, L.m);
  U.idiag = MALLOC(double, U.m);

  L.rowptr -= base;
  U.rowptr -= base;
  L.idiag -= base;
  U.idiag -= base;

  const int *perm = schedule.threadContToOrigPerm;

  int prefixSumWorkspace[2*(omp_get_max_threads() + 1)];
  int sum[2];

#pragma omp parallel
  {
    int tid = omp_get_thread_num();

    int iBegin, iEnd;
    getSimpleThreadPartition(&iBegin, &iEnd, A.m);
    iBegin += base; iEnd += base;

    // count # of nnz per row
    int sumPrivate[2] = { 0, 0 };
    for (int i = iBegin; i < iEnd; i++) {
      int row = perm[i];
      sumPrivate[0] += A.diagptr[row] - A.rowptr[row] + 1;
      sumPrivate[1] += A.rowptr[row + 1] - A.diagptr[row];
    }

    prefixSumMultiple(sumPrivate, sum, 2, prefixSumWorkspace);

#pragma omp master
    {
      L.rowptr[L.m + base] = sum[0] + base;
      U.rowptr[U.m + base] = sum[1] + base;

      L.colidx = MALLOC(int, sum[0]);
      U.colidx = MALLOC(int, sum[1]);

      L.values = MALLOC(double, sum[0]);
      U.values = MALLOC(double, sum[1]);

      L.colidx -= base;
      U.colidx -= base;

      L.values -= base;
      U.values -= base;
    }
#pragma omp barrier

    sumPrivate[0] += base;
    sumPrivate[1] += base;

    const int ntasks = schedule.ntasks;
    const short *nparents = schedule.nparentsForward;
    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    int nBegin, nEnd;
    getSimpleThreadPartition(&nBegin, &nEnd, ntasks);

    volatile int *taskFinished = schedule.taskFinished;
    int **parents = schedule.parentsForward;

    memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

    synk::Barrier::getInstance()->wait(tid);

    for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
      SPMP_LEVEL_SCHEDULE_WAIT;

      for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
        int row = perm[i] + base;

        for (int j = rowptr[row]; j < diagptr[row]; ++j) {
          int c = colidx[j];
          double tmp = lu[j] /= lu[diagptr[c]];

          int k1 = j + 1, k2 = diagptr[c] + 1;

          while (k1 < rowptr[row + 1] && k2 < rowptr[c + 1]) {
            if (colidx[k1] < colidx[k2]) ++k1;
            else if (colidx[k1] > colidx[k2]) ++k2;
            else {
              lu[k1] -= tmp*lu[k2];
              ++k1; ++k2;
            }
          }
        }
      } // for each row

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each level
  } // omp parallel

  L.rowptr += base;
  U.rowptr += base;

  L.idiag += base;
  U.idiag += base;

  L.colidx += base;
  U.colidx += base;

  L.values += base;
  U.values += base;
}*/

} // namespace SpMP
