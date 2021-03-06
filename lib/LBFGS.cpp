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
#include <cstring>
#include <algorithm>
#include <cmath> // fabs, exp, loglp

#include "CSR_Interface.h"

extern "C" {

void LBFGSComputeDirection(int k, int it, int n, const double *S, const double *Y, const double *dfk, double *r)
{
  double a[k];
  int mm = std::min(it - 1, k);
  int begin_idx = it + k - mm - 2;

  double yk = 1;

  double p[mm];

  if (mm > 0) {
    int j = (begin_idx + mm)%k;
    
    double sum = 0;
    double p_sum = 0;
    if (it > k) {
      double sum2 = 0;
#pragma omp parallel for reduction(+:sum,p_sum,sum2)
      for (int i = 0; i < n; ++i) {
        sum += S[j*n + i]*dfk[i];
        p_sum += S[j*n + i]*Y[j*n + i];
        sum2 += Y[j*n + i]*Y[j*n + i];
      }
      yk = p_sum/sum2;
    }
    else {
#pragma omp parallel for reduction(+:sum,p_sum)
      for (int i = 0; i < n; ++i) {
        sum += S[j*n + i]*dfk[i];
        p_sum += S[j*n + i]*Y[j*n + i];
      }
    }
    p[mm - 1] = 1/p_sum;
    a[mm - 1] = p[mm - 1]*sum;
  }

  for (int l = mm; l >= 2; --l) {
    int j = (begin_idx + l)%k;
    int j_minus_1 = (j + k - 1)%k;
    
    double sum = 0;
    double p_sum = 0;
    if (l == mm) {
#pragma omp parallel for reduction(+:sum,p_sum)
      for (int i = 0; i < n; ++i) {
        r[i] = dfk[i] - a[l - 1]*Y[j*n + i];
        sum += S[j_minus_1*n + i]*r[i];
        p_sum += S[j_minus_1*n + i]*Y[j_minus_1*n + i];
      }
    }
    else {
#pragma omp parallel for reduction(+:sum,p_sum)
      for (int i = 0; i < n; ++i) {
        r[i] -= a[l - 1]*Y[j*n + i];
        sum += S[j_minus_1*n + i]*r[i];
        p_sum += S[j_minus_1*n + i]*Y[j_minus_1*n + i];
      }
    }

    // a loop carried depdendency via a[:]
    p[l - 2] = 1/p_sum;
    a[l - 2] = p[l - 2]*sum;
  }

  if (mm > 0) {
    int j = (begin_idx + 1)%k;

    if (mm == 1) {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        r[i] = dfk[i] - a[0]*Y[j*n + i];
        r[i] = -yk*r[i];
      }
    }
    else {
#pragma omp parallel for
      for (int i = 0; i < n; ++i) {
        r[i] -= a[0]*Y[j*n + i];
        r[i] = -yk*r[i];
      }
    }
  }
  else {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      r[i] = -yk*dfk[i];
    }
  }

  // this loop can be parallelized by privatization of r
  for (int l = 1; l <= mm; ++l) {
    int j = (begin_idx + l)%k;

    double b = 0;
#pragma omp parallel for reduction(+:b)
    for (int i = 0; i < n; ++i) {
      b -= Y[j*n + i]*r[i];
    }
    b *= p[l - 1];
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      r[i] -= (a[l - 1] - b)*S[j*n + i];
    }
  }
}

// s = sum(log1p(exp(-abs(yXw)) - min(yXw,0)))
// fk = s/m + (lambda/2)*norm(w)^2
double LBFGSLossFunction1(
  int m, int n, const double *yXw, const double *w, double lambda)
{
  double s = 0;
#pragma omp parallel for reduction(+:s)
  for (int i = 0; i < m; ++i) {
    s += log1p(exp(-fabs(yXw[i]))) - std::min(yXw[i], 0.);
  }

  return s/m + (lambda/2)*dot(n, w, w);
}

// temp = y./(1 + exp(yXw))
void LBFGSLossFunction2(
  double *temp, int n, const double *y, const double *yXw)
{
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    temp[i] = y[i]/(1 + exp(yXw[i]));
  }
}

} // extern "C"
