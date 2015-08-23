#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "SpMP/CSR.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"
using namespace SpMP;

void forwardSolveRef(const CSR& A, double y[], const double b[]);

void backwardSolveRef(const CSR& A, double y[], const double b[]);

void forwardSolve(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm);
  
void backwardSolve(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm);

void forwardSolveWithReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule);
 
void backwardSolveWithReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule);

#endif
