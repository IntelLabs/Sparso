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
  const LevelSchedule& schedule,
  const int *perm)
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
  const LevelSchedule& schedule, const int *perm)
{
  if (0 == A.getBase()) {
    forwardSolve_<0>(A, y, b, schedule, perm);
  }
  else {
    assert(1 == A.getBase());
    forwardSolve_<1>(A, y, b, schedule, perm);
  }
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
template<int BASE = 0>
static void backwardSolve_(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm)
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
        int row = perm[i];
        double sum = b[row];
        for (int j = A.rowptr[row + 1] - 1 - BASE; j >= A.rowptr[row] - BASE; --j) {
          sum -= A.values[j]*y[A.colidx[j] - BASE];
        }
        y[row] += sum*A.idiag[i];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

void backwardSolve(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, const int *perm)
{
  if (0 == A.getBase()) {
    backwardSolve_<0>(A, y, b, schedule, perm);
  }
  else {
    assert(1 == A.getBase());
    backwardSolve_<1>(A, y, b, schedule, perm);
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
