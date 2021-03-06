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
#include "../BFSBipartite.h"
#include "SpMP/test/test.hpp"

using namespace SpMP;

int main(int argc, char **argv)
{
  CSR *A = new CSR(argv[1], 1);
  CSR *AT = A->transpose();

  int nnz = A->getNnz();
  double flops = 2*nnz;
  double bytes = (sizeof(double) + sizeof(int))*nnz + sizeof(double)*(A->m + A->n);

  printf("original bandwidth %d\n", A->getBandwidth());

  double *x = MALLOC(double, A->n);
  double *y = MALLOC(double, A->m);

  // allocate a large buffer to flush out cache
  bufToFlushLlc = (double *)_mm_malloc(LLC_CAPACITY, 64);

  const int REPEAT = 128;
  double times[REPEAT];

  for (int i = 0; i < REPEAT; ++i) {
    flushLlc();

    double t = omp_get_wtime();
    A->multiplyWithVector(y, x);
    times[i] = omp_get_wtime() - t;
  }

  printf("SpMV BW");
  printEfficiency(times, REPEAT, flops, bytes);

  for (int i = 0; i < REPEAT; ++i) {
    flushLlc();

    double t = omp_get_wtime();
    AT->multiplyWithVector(x, y);
    times[i] = omp_get_wtime() - t;
  }

  printf("SpMVT BW");
  printEfficiency(times, REPEAT, flops, bytes);

  int *rowPerm = new int[A->m];
  int *rowInversePerm = new int[A->m];
  int *colPerm = new int[A->n];
  int *colInversePerm = new int[A->n];

  bfsBipartite(*A, *AT, rowPerm, rowInversePerm, colPerm, colInversePerm);

  //for (int i = 0; i < A->n; ++i) {
    //colPerm[i] = i;
  //}

  CSR *APerm = A->permute(colPerm, rowInversePerm);
  CSR *ATPerm = A->transpose();

  printf("Permuted bandwidth %d\n", APerm->getBandwidth());

  for (int i = 0; i < REPEAT; ++i) {
    flushLlc();

    double t = omp_get_wtime();
    APerm->multiplyWithVector(y, x);
    times[i] = omp_get_wtime() - t;
  }

  printf("SpMV BW");
  printEfficiency(times, REPEAT, flops, bytes);

  for (int i = 0; i < REPEAT; ++i) {
    flushLlc();

    double t = omp_get_wtime();
    ATPerm->multiplyWithVector(x, y);
    times[i] = omp_get_wtime() - t;
  }

  printf("SpMV BW");
  printEfficiency(times, REPEAT, flops, bytes);

  return 0;
}
