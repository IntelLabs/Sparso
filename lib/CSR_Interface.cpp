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
#include "SpMP/CSR.hpp"
#include "CSR_Interface.h"
#include "SpGEMM.hpp"
#include <assert.h>

#define PERF_TUNE

#ifdef PERF_TUNE
#include <sys/time.h>
#include <stdio.h>
#include <omp.h>
#endif

using namespace SpMP;

extern "C" {

CSR_Handle *CSR_Create(int numRows, int numCols, int *i, int *j, double *v)
{
  CSR *ret = new CSR;

  ret->m = numRows;
  ret->n = numCols;
  ret->rowptr = i;
  ret->colidx = j;
  ret->values = v;
  
  return (CSR_Handle *)ret;
}

int CSR_GetNumRows(CSR_Handle *A)
{
  return ((CSR *)A)->m;
}

int CSR_GetNumCols(CSR_Handle *A)
{
  return ((CSR *)A)->n;
}

int CSR_GetNumNonZeros(CSR_Handle *A)
{
  return ((CSR *)A)->rowptr[((CSR *)A)->m];
}

int *CSR_GetRowPtr(CSR_Handle *A)
{
  return ((CSR *)A)->rowptr;
}

int *CSR_GetColIdx(CSR_Handle *A)
{
  return ((CSR *)A)->colidx;
}

double *CSR_GetValues(CSR_Handle *A)
{
  return ((CSR *)A)->values;
}

CSR_Handle *CSR_ADBInspect(const CSR_Handle *A, const CSR_Handle *B)
{
    return (CSR_Handle *)SpMP::inspectADB((CSR*)A, (CSR*)B);
}

void CSR_ADB(CSR_Handle *C, const CSR_Handle *A, const CSR_Handle *B, const double *d)
{
    adb((CSR *)C, (CSR *)A, (CSR *)B, (double *)d);
}

void CSR_GetRCMPermutation(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm);
}

void CSR_GetRCMPermutationWithoutPseudoDiameterSourceSelection(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm, false);
}

void CSR_GetBFSPermutation(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->getBFSPermutation(perm, inversePerm);
}

void CSR_Permute(const CSR_Handle *A, CSR_Handle *out, const int *columnPerm, const int *rowInversePerm)
{
  ((CSR *)A)->permuteRowptr((CSR *)out, rowInversePerm);
  ((CSR *)A)->permuteMain((CSR *)out, columnPerm, rowInversePerm);
}

int CSR_GetBandwidth(CSR_Handle *A)
{
  return ((CSR *)A)->getBandwidth();
}

void CSR_Destroy(CSR_Handle *A)
{
  delete (CSR *)A;
}

#ifdef PERF_TUNE
static double* stats = NULL;
void CSR_Statistics(double *_stats)
{
    stats = _stats;
}
#endif

void CSR_Make0BasedIndexing(CSR_Handle *A)
{
  ((CSR *)A)->make0BasedIndexing();
}

void CSR_Make1BasedIndexing(CSR_Handle *A)
{
  ((CSR *)A)->make1BasedIndexing();
}

} // extern "C"
