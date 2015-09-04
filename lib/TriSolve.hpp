#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "SpMP/CSR.hpp"
#include "SpMP/LevelSchedule.hpp"

using namespace SpMP;

namespace SpMP
{

void forwardSolveRef(CSR& A, double y[], const double b[]);

void backwardSolveRef(CSR& A, double y[], const double b[]);

void forwardSolve(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm);
  
void backwardSolve(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm);

void forwardSolveWithReorderedMatrix(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule);
 
void backwardSolveWithReorderedMatrix(
  CSR& A, double y[], const double b[],
  const LevelSchedule& schedule);

} // namespace SpMP

#endif
