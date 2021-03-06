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
#include <algorithm>
#include <cassert>
#include <cstring>
#include <climits>
#include <cfloat>
#include <omp.h>
#include "TriSolve.hpp"
#include "CSR_Interface.h"
#include "SpMP/CSR.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"

using namespace std;

namespace SpMP
{

/**
 * Reference sequential sparse triangular solver
 */
template<int BASE = 0>
static void forwardSolveRef_(CSR& A, double y[], const double b[])
{
  A.computeInverseDiag();

  for (int i = 0; i < A.m; ++i) {
    double sum = b[i];
    for (int j = A.rowptr[i] - BASE; j < A.rowptr[i + 1] - BASE; ++j) {
      sum -= A.values[j]*y[A.colidx[j] - BASE];
    }
    y[i] += sum*A.idiag[i];
  } // for each row
}

void forwardSolveRef(CSR& A, double y[], const double b[])
{
  if (0 == A.getBase()) {
    forwardSolveRef_<0>(A, y, b);
  }
  else {
    assert(1 == A.getBase());
    forwardSolveRef_<1>(A, y, b);
  }
}

template<int BASE = 0>
static void backwardSolveRef_(CSR& A, double y[], const double b[])
{
  A.computeInverseDiag();

  for (int i = A.m - 1; i >= 0; --i) {
    double sum = b[i];
    for (int j = A.rowptr[i] - BASE; j < A.rowptr[i + 1] - BASE; ++j) {
      sum -= A.values[j]*y[A.colidx[j] - BASE];
    }
    y[i] += sum*A.idiag[i];
  } // for each row
}

void backwardSolveRef(CSR& A, double y[], const double b[])
{
  if (0 == A.getBase()) {
    backwardSolveRef_<0>(A, y, b);
  }
  else {
    assert(1 == A.getBase());
    backwardSolveRef_<1>(A, y, b);
  }
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void forwardSolve_(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  A.computeInverseDiag();
  const int *perm = schedule.threadContToOrigPerm;

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    const int ntasks = schedule.ntasks;
    const short *nparents = schedule.nparentsForward;
    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    int nPerThread = (ntasks + nthreads - 1)/nthreads;
    int nBegin = min(nPerThread*tid, ntasks);
    int nEnd = min(nBegin + nPerThread, ntasks);

    volatile int *taskFinished = schedule.taskFinished;
    int **parents = schedule.parentsForward;

    memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

    synk::Barrier::getInstance()->wait(tid);

    for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
      SPMP_LEVEL_SCHEDULE_WAIT;

      for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
        int row = perm[i];
        double sum = b[row];
        for (int j = A.rowptr[row] - BASE; j < A.rowptr[row + 1] - BASE; ++j) {
          sum -= A.values[j]*y[A.colidx[j] - BASE];
        }
        y[row] += sum*A.idiag[row];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void forwardSolve(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  if (0 == A.getBase()) {
    forwardSolve_<0>(A, y, b, schedule);
  }
  else {
    assert(1 == A.getBase());
    forwardSolve_<1>(A, y, b, schedule);
  }
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void backwardSolve_(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  A.computeInverseDiag();
  const int *perm = schedule.threadContToOrigPerm;

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    const int ntasks = schedule.ntasks;
    const short *nparents = schedule.nparentsBackward;
    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    int nPerThread = (ntasks + nthreads - 1)/nthreads;
    int nBegin = min(nPerThread*tid, ntasks);
    int nEnd = min(nBegin + nPerThread, ntasks);

    volatile int *taskFinished = schedule.taskFinished;
    int **parents = schedule.parentsBackward;

    memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

    synk::Barrier::getInstance()->wait(tid);

    for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task) {
      SPMP_LEVEL_SCHEDULE_WAIT;

      for (int i = taskBoundaries[task + 1] - 1; i >= taskBoundaries[task]; --i) {
        int row = perm[i];
        double sum = b[row];
        for (int j = A.rowptr[row + 1] - 1 - BASE; j >= A.rowptr[row] - BASE; --j) {
          sum -= A.values[j]*y[A.colidx[j] - BASE];
        }
        y[row] += sum*A.idiag[row];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void backwardSolve(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  if (0 == A.getBase()) {
    backwardSolve_<0>(A, y, b, schedule);
  }
  else {
    assert(1 == A.getBase());
    backwardSolve_<1>(A, y, b, schedule);
  }
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void forwardSolveWithReorderedMatrix_(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  A.computeInverseDiag();

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    const int ntasks = schedule.ntasks;
    const short *nparents = schedule.nparentsForward;
    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    int nPerThread = (ntasks + nthreads - 1)/nthreads;
    int nBegin = min(nPerThread*tid, ntasks);
    int nEnd = min(nBegin + nPerThread, ntasks);

    volatile int *taskFinished = schedule.taskFinished;
    int **parents = schedule.parentsForward;

    memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

    synk::Barrier::getInstance()->wait(tid);

    for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
      SPMP_LEVEL_SCHEDULE_WAIT;

      for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
        double sum = b[i];
        for (int j = A.rowptr[i] - BASE; j < A.rowptr[i + 1] - BASE; ++j) {
          sum -= A.values[j]*y[A.colidx[j] - BASE];
        }
        y[i] += sum*A.idiag[i];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void forwardSolveWithReorderedMatrix(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  if (0 == A.getBase()) {
    forwardSolveWithReorderedMatrix_<0>(A, y, b, schedule);
  }
  else {
    assert(1 == A.getBase());
    forwardSolveWithReorderedMatrix_<1>(A, y, b, schedule);
  }
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void backwardSolveWithReorderedMatrix_(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  A.computeInverseDiag();

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    const int ntasks = schedule.ntasks;
    const short *nparents = schedule.nparentsBackward;
    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    int nPerThread = (ntasks + nthreads - 1)/nthreads;
    int nBegin = min(nPerThread*tid, ntasks);
    int nEnd = min(nBegin + nPerThread, ntasks);

    volatile int *taskFinished = schedule.taskFinished;
    int **parents = schedule.parentsBackward;

    memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

    synk::Barrier::getInstance()->wait(tid);

    for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task) {
      SPMP_LEVEL_SCHEDULE_WAIT;

      for (int i = taskBoundaries[task + 1] - 1; i >= taskBoundaries[task]; --i) {
        double sum = b[i];
        for (int j = A.rowptr[i + 1] - 1 - BASE; j >= A.rowptr[i] - BASE; --j) {
          sum -= A.values[j]*y[A.colidx[j] - BASE];
        }
        y[i] += sum*A.idiag[i];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void backwardSolveWithReorderedMatrix(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  if (0 == A.getBase()) {
    backwardSolveWithReorderedMatrix_<0>(A, y, b, schedule);
  }
  else {
    assert(1 == A.getBase());
    backwardSolveWithReorderedMatrix_<1>(A, y, b, schedule);
  }
}

} // namespace SpMP
