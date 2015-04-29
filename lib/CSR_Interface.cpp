#include "CSR.hpp"
#include "CSR_Interface.h"
#include <assert.h>

//#define PERF_TUNE

#ifdef PERF_TUNE
#include <sys/time.h>
#include <stdio.h>
#include <omp.h>
#endif

extern "C" {

CSR_Handle *CSR_Create(int numRows, int numCols, int *i, int *j, double *v, int base)
{
  return (CSR_Handle *)(new CSR(numRows, numCols, i, j, v, base));
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
  return ((CSR *)A)->rowPtr[((CSR *)A)->m];
}

void CSR_MultiplyWithVector(const CSR_Handle *A, double *y, const double *x)
{
  ((CSR *)A)->multiplyWithVector(y, x);
}

#ifdef USE_BOOST
void CSR_BoostGetRCMPemutation(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->boostGetRCMPermutation(perm, inversePerm);
}

void CSR_BoostGetRCMPemutationWithSource(const CSR_Handle *A, int *perm, int *inversePerm, int source)
{
  ((CSR *)A)->boostGetRCMPermutation(perm, inversePerm, source);
}
#endif

void CSR_GetRCMPemutation(const CSR_Handle *A, int *perm, int *inversePerm)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm);
}

void CSR_GetRCMPemutationWithSource(const CSR_Handle *A, int *perm, int *inversePerm, int source)
{
  ((CSR *)A)->getRCMPermutation(perm, inversePerm, source);
}

void CSR_Permute(const CSR_Handle *A, CSR_Handle *out, const int *columnPerm, const int *rowInversePerm)
{
  ((CSR *)A)->permute((CSR *)out, columnPerm, rowInversePerm);
}

int CSR_GetBandwidth(CSR_Handle *A)
{
  return ((CSR *)A)->getBandwidth();
}

void CSR_PrintInDense(CSR_Handle *A)
{
  ((CSR *)A)->printInDense();
}

void CSR_PrintSomeValues(int numRows, int numCols, int *i, int *j, double *v, int distance, bool is_1_based)
{
  CSR *A = new CSR(numRows, numCols, i, j, v);
  A->printSomeValues(distance, is_1_based);
}

void CSR_Destroy(CSR_Handle *A)
{
  delete (CSR *)A;
}

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
#ifdef PERF_TUNE
    double t1 = omp_get_wtime();
#endif

    // The original and the result array space must be different
    assert(i != i1);
    assert(j != j1);    
    assert(v != v1);
    
    CSR *A = new CSR(numRows, numCols, i, j, v, oneBasedOutput ? 1 : 0);
    if (oneBasedInput) {
        // The library assumes arrays are all 0-based indexing
        A->make0BasedIndexing();
    }

#ifdef PERF_TUNE
    double t2 = omp_get_wtime();
#endif

    if (getPermutation) {
        A->getRCMPermutation(perm, inversePerm);
    }

#ifdef PERF_TUNE
    double t3 = omp_get_wtime();
#endif

    CSR *newA = new CSR(numRows, numCols, i1, j1, v1);
    A->permute(newA, perm, inversePerm);
    
#ifdef PERF_TUNE
    double t4 = omp_get_wtime();
#endif

    if (oneBasedOutput) {
        newA->make1BasedIndexing();
    }
    if (oneBasedInput) {
        // Recover the inpt from 0-based to 1-based indexing
        A->make1BasedIndexing();
    }
    delete newA;
    delete A;

#ifdef PERF_TUNE
    double t5 = omp_get_wtime();

    printf("CSR_ReorderMatrix total: %f sec\n", t5 - t1);
    printf("\tmake 0 based: %f sec\n", t2 - t1);
    printf("\tCSR_GetRCMPemutation: %f sec\n", t3 - t2);
    printf("\tCSR_Permute: %f sec\n", t4 - t3);
    printf("\tmake 1 based: %f sec\n", t5 - t4);
    fflush(stdout);
#endif
}

}
