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
void forwardSolveRef(const CSR& A, double y[], const double b[])
{
  for (int i = 0; i < A.m; ++i) {
    double sum = b[i];
    for (int j = A.rowptr[i]; j < A.rowptr[i + 1]; ++j) {
      sum -= A.values[j]*y[A.colidx[j]];
    }
    y[i] = sum*A.idiag[i];
  } // for each row
}

void backwardSolveRef(const CSR& A, double y[], const double b[])
{
  for (int i = A.m - 1; i >= 0; --i) {
    double sum = b[i];
    for (int j = A.rowptr[i]; j < A.rowptr[i + 1]; ++j) {
      sum -= A.values[j]*y[A.colidx[j]];
    }
    y[i] = sum;
  } // for each row
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization
 */
void forwardSolveWithBarrier(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm)
{
#pragma omp parallel
  {
    int tid = omp_get_thread_num();

    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
      for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
        int row = perm[i];
        double sum = b[row];
        for (int j = A.rowptr[row]; j < A.rowptr[row + 1]; ++j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[row] = sum*A.idiag[row];
      } // for each row

      synk::Barrier::getInstance()->wait(tid);
    } // for each level
  } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization
 */
void backwardSolveWithBarrier(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm)
{
#pragma omp parallel
  {
    int tid = omp_get_thread_num();

    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task) {
      for (int i = taskBoundaries[task + 1] - 1; i >= taskBoundaries[task]; --i) {
        int row = perm[i];
        double sum = b[row];
        for (int j = A.rowptr[row + 1] - 1; j >= A.rowptr[row]; --j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[row] = sum;
      } // for each row
      synk::Barrier::getInstance()->wait(tid);
    } // for each level
  } // omp parallel
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void forwardSolve(
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
        for (int j = A.rowptr[row]; j < A.rowptr[row + 1]; ++j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[row] = sum*A.idiag[row];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void backwardSolve(
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
        for (int j = A.rowptr[row + 1] - 1; j >= A.rowptr[row]; --j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[row] = sum;
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization. Matrix is reordered.
 */
void forwardSolveWithBarrierAndReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
#pragma omp parallel
  {
    int tid = omp_get_thread_num();

    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
      for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
        double sum = b[i];
        for (int j = A.rowptr[i]; j < A.rowptr[i + 1]; ++j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[i] = sum*A.idiag[i];
      } // for each row
      synk::Barrier::getInstance()->wait(tid);
    } // for each level
  } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization. Matrix is reordered.
 */
void backwardSolveWithBarrierAndReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule)
{
#pragma omp parallel
  {
    int tid = omp_get_thread_num();

    const vector<int>& threadBoundaries = schedule.threadBoundaries;
    const vector<int>& taskBoundaries = schedule.taskBoundaries;

    for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task) {
      for (int i = taskBoundaries[task + 1] - 1; i >= taskBoundaries[task]; --i) {
        double sum = b[i];
        for (int j = A.rowptr[i + 1] - 1; j >= A.rowptr[i]; --j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[i] = sum;
      } // for each row
      synk::Barrier::getInstance()->wait(tid);
    } // for each level
  } // omp parallel
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void forwardSolveWithReorderedMatrix(
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
        for (int j = A.rowptr[i]; j < A.rowptr[i + 1]; ++j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[i] = sum*A.idiag[i];
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void backwardSolveWithReorderedMatrix(
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
        for (int j = A.rowptr[i + 1] - 1; j >= A.rowptr[i]; --j) {
          sum -= A.values[j]*y[A.colidx[j]];
        }
        y[i] = sum;
      }

      SPMP_LEVEL_SCHEDULE_NOTIFY;
    } // for each task
  } // omp parallel
}
