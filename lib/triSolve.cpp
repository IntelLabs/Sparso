#include <algorithm>
#include <cassert>
#include <cstring>
#include <climits>
#include <cfloat>
#include <omp.h>
#include "triSolve.hpp"
#include "CSR_Interface.h"
#include "SpMP/CSR.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"

using namespace std;
using namespace SpMP;

/**
 * Reference sequential sparse triangular solver
 */
template<int BASE = 0>
static void forwardSolveRef_(const CSR& A, double y[], const double b[])
{
  for (int i = 0; i < A.m; ++i) {
    double sum = b[i];
    for (int j = A.rowptr[i] - BASE; j < A.rowptr[i + 1] - BASE; ++j) {
      sum -= A.values[j]*y[A.colidx[j] - BASE];
    }
    y[i] = sum*A.idiag[i];
  } // for each row
}

void forwardSolveRef(const CSR& A, double y[], const double b[])
{
  if (0 == A.base) {
    forwardSolveRef_<0>(A, y, b);
  }
  else {
    assert(1 == A.base);
    forwardSolveRef_<1>(A, y, b);
  }
}

template<int BASE = 0>
static void backwardSolveRef_(const CSR& A, double y[], const double b[])
{
  for (int i = A.m - 1; i >= 0; --i) {
    double sum = b[i];
    for (int j = A.rowptr[i] - BASE; j < A.rowptr[i + 1] - BASE; ++j) {
      sum -= A.values[j]*y[A.colidx[j] - BASE];
    }
    y[i] = sum;
  } // for each row
}

void backwardSolveRef(const CSR& A, double y[], const double b[])
{
  if (0 == A.base) {
    backwardSolveRef_<0>(A, y, b);
  }
  else {
    assert(1 == A.base);
    backwardSolveRef_<1>(A, y, b);
  }
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void forwardSolve_(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm)
{
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
        y[row] = sum*A.idiag[row];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void forwardSolve(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, const int *perm)
{
  if (0 == A.base) {
    forwardSolve_<0>(A, y, b, schedule, perm);
  }
  else {
    assert(1 == A.base);
    forwardSolve_<1>(A, y, b, schedule, perm);
  }
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void backwardSolve_(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm)
{
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
        y[row] = sum;
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void backwardSolve(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, const int *perm)
{
  if (0 == A.base) {
    backwardSolve_<0>(A, y, b, schedule, perm);
  }
  else {
    assert(1 == A.base);
    backwardSolve_<1>(A, y, b, schedule, perm);
  }
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void forwardSolveWithReorderedMatrix_(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
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
        y[i] = sum*A.idiag[i];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void forwardSolveWithReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  if (0 == A.base) {
    forwardSolveWithReorderedMatrix_<0>(A, y, b, schedule);
  }
  else {
    assert(1 == A.base);
    forwardSolveWithReorderedMatrix_<1>(A, y, b, schedule);
  }
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void backwardSolveWithReorderedMatrix_(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
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
        y[i] = sum;
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void backwardSolveWithReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
  if (0 == A.base) {
    backwardSolveWithReorderedMatrix_<0>(A, y, b, schedule);
  }
  else {
    assert(1 == A.base);
    backwardSolveWithReorderedMatrix_<1>(A, y, b, schedule);
  }
}
