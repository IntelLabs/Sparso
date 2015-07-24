#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "SpMP/CSR.hpp"
#include "SpMP/LevelSchedule.hpp"
#include "SpMP/synk/barrier.hpp"
using namespace SpMP;

void forwardSolveRef(const CSR& A, double y[], const double b[]);

void backwardSolveRef(const CSR& A, double y[], const double b[]);

void forwardSolveWithBarrier(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm, synk::Barrier *bar);

void backwardSolveWithBarrier(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm, synk::Barrier *bar);

void forwardSolve(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm, synk::Barrier *bar);
  
void backwardSolve(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule,
  const int *perm, synk::Barrier *bar);

void forwardSolveWithBarrierAndReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, synk::Barrier *bar);

void backwardSolveWithBarrierAndReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, synk::Barrier *bar);
 
void forwardSolveWithReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, synk::Barrier *bar);
 
void backwardSolveWithReorderedMatrix(
  const CSR& A, double y[], const double b[],
  const LevelSchedule& schedule, synk::Barrier *bar);

#endif
