#include "SpMP/CSR.hpp"
#include "CSR_Interface.h"
#include <assert.h>

#define PERF_TUNE

#ifdef PERF_TUNE
#include <sys/time.h>
#include <stdio.h>
#include <omp.h>
#endif

using namespace SpMP;

extern "C" {

CSR_Handle *CSR_Create(int numRows, int numCols, int *i, int *j, double *v, int base)
{
  CSR *ret = new CSR;

  ret->m = numRows;
  ret->n = numCols;
  ret->rowptr = i;
  ret->colidx = j;
  ret->values = v;
  ret->base = base;
  
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

// The first few parameters (numRows to v) represent the source matrix to be reordered. The next few
// parameters (i1~v1) are the spaces that have been allocated to store the results. 
// Perm and inversePerm are the spaces that have been allocated for permutation and inverse permutation
// info; when getPermutation is true, this function computes and stores the info into them; otherwise,
// they already contain the info, and this function just uses it.
// oneBasedInput: true if the i and j are 1-based indexing.
// oneBasedOutput: true if i1 and j1 should be 1-based indexing  
void CSR_ReorderMatrix(int numRows, int numCols, int *i, int *j, double *v, int *i1, int *j1, double *v1, 
                 int *perm, int *inversePerm, bool getPermutation, bool oneBasedInput, bool oneBasedOutput)
{
    double t1, t2, t3, t4, t5;
#ifdef PERF_TUNE
    t1 = omp_get_wtime();
#endif

    // The original and the result array space must be different
    assert(i != i1);
    assert(j != j1);    
    assert(v != v1);
    
    CSR *A = new CSR(numRows, numCols, i, j, v, oneBasedInput ? 1 : 0);

#ifdef PERF_TUNE
    int orig_bw = A->getBandwidth();
    t2 = omp_get_wtime();
#endif

    if (getPermutation) {
        A->getBFSPermutation(perm, inversePerm);
    }

#ifdef PERF_TUNE
    t3 = omp_get_wtime();
#endif

    CSR *newA = new CSR(numRows, numCols, i1, j1, v1, oneBasedInput ? 1 : 0);
    A->permuteRowptr(newA, inversePerm);
    A->permuteMain(newA, perm, inversePerm);
    
#ifdef PERF_TUNE
    int rcm_bw = newA->getBandwidth();
    t4 = omp_get_wtime();
#endif

    if (oneBasedOutput) {
        newA->make1BasedIndexing();
    }
    else {
        newA->make0BasedIndexing();
    }
    delete newA;
    delete A;

#ifdef PERF_TUNE
        t5 = omp_get_wtime();
        double bytes = (double)i[numRows]*12;
        if (stats) {
          stats[0] += t5 - t1;
          stats[1] += t3 - t2;
          stats[2] += t4 - t3;
        }
    
#if 0
        printf("CSR_ReorderMatrix total: %f sec\n", t5 - t1);
        printf("\tCSR_GetRCMPermutation: %f sec (%f GB/s)\n", t3 - t2, bytes/(t3 - t2)/1e9);
        printf("\tCSR_Permute: %f sec (%f GB/s)\n", t4 - t3, bytes/(t4 - t3)/1e9);
        printf("\tBW changed: %d -> %d\n", orig_bw, rcm_bw);
        fflush(stdout);
#endif

#endif
}

void CSR_Make0BasedIndexing(CSR_Handle *A)
{
  ((CSR *)A)->make0BasedIndexing();
}

void CSR_Make1BasedIndexing(CSR_Handle *A)
{
  ((CSR *)A)->make1BasedIndexing();
}

} // extern "C"
