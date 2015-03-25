#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <algorithm>
#include <omp.h>

#include "CSR_Interface.h"
#include "mm_io.h"

#define T double

static void sort(int *col_idx, T *a, int start, int end)
{
  int i, j, it;
  double dt;

  for (i=end-1; i>start; i--)
    for(j=start; j<i; j++)
      if (col_idx[j] > col_idx[j+1]){

  if (a){
    dt=a[j]; 
    a[j]=a[j+1]; 
    a[j+1]=dt;
        }
  it=col_idx[j]; 
  col_idx[j]=col_idx[j+1]; 
  col_idx[j+1]=it;
    
      }
}

/* converts COO format (1-based index) to CSR format (0-based index), not in-place,*/
static void coo2csr(int n, int nz, T *a, int *i_idx, int *j_idx,
       T *csr_a, int *col_idx, int *row_start)
{
  int i, l;

  for (i=0; i<=n; i++) row_start[i] = 0;

  /* determine row lengths */
  for (i=0; i<nz; i++) row_start[i_idx[i]]++;


  for (i=0; i<n; i++) row_start[i+1] += row_start[i];


  /* go through the structure  once more. Fill in output matrix. */
  for (l=0; l<nz; l++){
    i = row_start[i_idx[l] - 1];
    csr_a[i] = a[l];
    col_idx[i] = j_idx[l] - 1;
    row_start[i_idx[l] - 1]++;
  }

  /* shift back row_start */
  for (i=n; i>0; i--) row_start[i] = row_start[i-1];

  row_start[0] = 0;

  for (i=0; i<n; i++){
    sort (col_idx, csr_a, row_start[i], row_start[i+1]);
  }

}

void set_1based_ind(int *rowptr, int *colidx, int n, int nnz)
{
  int i;
  for(i=0; i <= n; i++)
     rowptr[i]++;
  for(i=0; i < nnz; i++)
     colidx[i]++;
}

void set_0based_ind(int *rowptr, int *colidx, int n, int nnz)
{
  int i;
  for(i=0; i <= n; i++)
     rowptr[i]--;
  for(i=0; i < nnz; i++)
     colidx[i]--;
}


static double randfp(double low, double high){
    double t = (double)rand() / (double)RAND_MAX;
    return (1.0 - t) * low + t * high;
}

#define ZERO_BASED_INDEX
#define RAND01() randfp(0.0, 1.0)
void load_matrix_market (char *file, T **a, int **aj, int **ai, int *am, int *an, int *annz)
{
    FILE *fp=fopen(file, "r");
    assert(fp);
    MM_typecode matcode;
    int m;
    int n;
    int nnz;
    int x;
    int y;
    double value;
    int count;
    int pattern;
    bool symm = false;
    int i;
    int *colidx;
    int *rowidx;
    T *values;
    int lines;

    if (mm_read_banner (fp, &matcode) != 0)
    {
        printf ("Error: could not process Matrix Market banner.\n");
        exit(1);
    }

    if ( !mm_is_valid (matcode) &&
         (mm_is_array (matcode) ||
          mm_is_dense (matcode)) )
    {
        printf ("Error: only support sparse and real matrices.\n");
        exit(1);
    }
    
    if (mm_read_mtx_crd_size(fp, &m, &n, &nnz) !=0)
    {
        printf ("Error: could not read matrix size.\n");
        exit(1);
    
    }

    if (mm_is_symmetric (matcode) == 1)
    {
        symm = true;
        count = 2*nnz;
    }
    else
    {
        count = nnz;
    }

    values = (T *)malloc (count * sizeof(T));
    colidx = (int *)malloc (count * sizeof(int));
    rowidx = (int *)malloc (count * sizeof(int));
    assert (values != NULL);
    assert (colidx != NULL);
    assert (rowidx != NULL);
    
    count = 0;
    lines = 0;
    int trc=0;
    pattern = mm_is_pattern (matcode);   
    int x_o=1, y_o;
    bool zero_based_ind = false;
    while (mm_read_mtx_crd_entry (fp, &x, &y, &value, NULL, matcode) == 0)
    {
        if (0 == x || 0 == y) zero_based_ind = true;

#ifdef IGNORE_ZERO_WHEN_READ
        if (0 == value) continue; // not sure if this is the right thing to do
#endif
        rowidx[count] = x;
        colidx[count] = y;

        if (x > m || y > n)
        {
            printf ("Error: (%d %d) coordinate is out of range.\n", x, y);
            exit(1);
        }
        if (pattern == 1)
        {
            values[count] = RAND01();
        }
        else
        {
            values[count] = value;
        }

        #if defined(DIAG_DOMINANT)
        if(x==y) values[count]*=1.14;
        #endif

        count++;
        lines++;

        if (symm && x != y)
        {
            rowidx[count] = y;
            colidx[count] = x;
            if (pattern == 1)
            {
                values[count] = RAND01();
            }
            else
            {
                values[count] = value;
            }
            count++;
            trc++;
        }

    }
    assert (lines <= nnz);
    nnz = lines;

    // convert to crs
    *a = (T *)_mm_malloc(sizeof(T)*nnz, 64);
    *aj = (int *)_mm_malloc(sizeof(int)*nnz, 64);
    *ai = (int *)_mm_malloc(sizeof(int)*m, 64);
    *am = m;
    *an = n;
    *annz = nnz;

    assert(!zero_based_ind);
    coo2csr(m, nnz, values, rowidx, colidx, *a, *aj, *ai);
#ifndef ZERO_BASED_INDEX
    set_1based_ind(*ai, *aj, m, nnz);
#endif

    free (colidx);
    free (rowidx);
    free (values);
    fclose(fp);

    const char *tag=(symm?"symmetric":"general");
    printf("%s:::%s m=%d n=%d nnz=%d\n", file, tag, m, n, nnz);
}
bool isPerm(int *perm, int n)
{
  int *temp = new int[n];
  memcpy(temp, perm, sizeof(int)*n);
  std::sort(temp, temp + n);
  int *last = std::unique(temp, temp + n);
  if (last != temp + n) {
    memcpy(temp, perm, sizeof(int)*n);
    std::sort(temp, temp + n);

    for (int i = 0; i < n; ++i) {
      if (temp[i] == i - 1) {
        printf("%d duplicated\n", i);
        assert(false);
        return false;
      }
      else if (temp[i] != i) {
        printf("%d missed\n", i);
        assert(false);
        return false;
      }
    }
  }
  delete[] temp;
  return true;
}

/* Expected output

original banwdith: 8
RCM permutation
0 8 5 7 3 6 4 2 1 9
RCM permuted bandwidth: 4

original banwdith: 6
RCM permutation
3 6 0 4 2 1 7 5
RCM permuted bandwidth: 2
 */
int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: CSR_test matrix_market_file\n");
    return -1;
  }

  {
    int m, n, nnz;
    double *a;
    int *aj, *ai;
    load_matrix_market(argv[1], &a, &aj, &ai, &m, &n, &nnz);
    printf("m = %d, n = %d, nnz = %d\n", m, n, nnz);

    CSR_Handle *A = CSR_Create(m, n, ai, aj, a);
    printf("original bandwidth: %d\n", CSR_GetBandwidth(A));

    double *x = (double *)malloc(sizeof(double)*n);
    double *y = (double *)malloc(sizeof(double)*m);

    const int REPEAT = 128;

    double t = -omp_get_wtime();
    for (int i = 0; i < REPEAT; ++i) {
      CSR_MultiplyWithVector(A, y, x);
    }
    t += omp_get_wtime();

    printf("SpMV BW = %g GB/s\n", ((double)nnz*12 + (m + n)*8)/(t/REPEAT)/1e9);

    printf("RCM permutation\n");
    int *perm = (int *)malloc(sizeof(int)*m);
    int *inversePerm = (int *)malloc(sizeof(int)*m);
    int source = 1;
    CSR_GetRCMPemutationWithSource(A, perm, inversePerm, source);
    isPerm(perm, m);
    isPerm(inversePerm, m);

    double *a2 = (double *)malloc(sizeof(double)*nnz);
    int *aj2 = (int *)malloc(sizeof(int)*nnz);
    int *ai2 = (int *)malloc(sizeof(int)*(m + 1));
    CSR_Handle *A2 = CSR_Create(m, n, ai2, aj2, a2);

    CSR_Permute(A, A2, perm, inversePerm);
    printf("RCM permuted bandwidth with source %d: %d\n\n", source, CSR_GetBandwidth(A2));

    CSR_GetRCMPemutation(A, perm, inversePerm);
    isPerm(perm, m);
    isPerm(inversePerm, m);
    
    CSR_Permute(A, A2, perm, inversePerm);
    printf("RCM permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

    CSR_GetRCMPemutationNew(A, perm, inversePerm, source);
    isPerm(perm, m);
    isPerm(inversePerm, m);
    
    CSR_Permute(A, A2, perm, inversePerm);
    printf("My RCM permuted bandwidth with source %d: %d\n\n", source, CSR_GetBandwidth(A2));

    t = -omp_get_wtime();
    for (int i = 0; i < REPEAT; ++i) {
      CSR_MultiplyWithVector(A2, y, x);
    }
    t += omp_get_wtime();

    printf("SpMV BW = %g GB/s\n", ((double)nnz*12 + (m + n)*8)/(t/REPEAT)/1e9);

    CSR_Destroy(A2);
    CSR_Destroy(A);
  }

#if 0
  {
    // from http://www.boost.org/doc/libs/1_37_0/libs/graph/example/cuthill_mckee_ordering.cpp
    int colIdx[] =    { 3, 5, 2, 4, 6, 9, 1, 3, 4, 0, 2, 5, 8, 1, 2, 6, 0, 3, 6, 7, 1, 4, 5, 7, 5, 6, 3, 1 };
    int rowPtr[] =    { 0,    2,          6,       9,          13,      16,         20,         24,   26,27, 28 };
    double values[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    int m = sizeof(rowPtr)/sizeof(rowPtr[0]) - 1;
    assert(10 == m);
    int nnz = rowPtr[m];
    assert(nnz == sizeof(colIdx)/sizeof(colIdx[0]));
    CSR_Handle *A = CSR_Create(m, m, rowPtr, colIdx, values);
    CSR_PrintInDense(A);
    printf("original banwdith: %d\n", CSR_GetBandwidth(A));

    printf("RCM permutation\n");
    int perm[m], inversePerm[m];
    CSR_GetRCMPemutation(A, perm, inversePerm);
    for (int i = 0; i < m; ++i) {
      printf("%d ", inversePerm[i]);
    }
    printf("\n");

    int rowPtr2[m + 1];
    int colIdx2[nnz];
    double values2[nnz];
    CSR_Handle *A2 = CSR_Create(m, m, rowPtr2, colIdx2, values2);

    CSR_Permute(A, A2, perm, inversePerm);
    CSR_PrintInDense(A2);
    printf("RCM permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

    CSR_Destroy(A2);
    CSR_Destroy(A);
  }

  {
    // from http://ciprian-zavoianu.blogspot.com/2009/01/project-bandwidth-reduction.html
    int colIdx[] =    { 0, 4, 1, 2, 5, 7, 1, 2, 4, 3, 6, 0, 2, 4, 1, 5, 7, 3, 6, 1, 5, 7 };
    int rowPtr[] =    { 0,    2,          6,       9,    11,      14,      17,   19,     22 };
    double values[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    int m = sizeof(rowPtr)/sizeof(rowPtr[0]) - 1;
    int nnz = rowPtr[m];
    assert(nnz == sizeof(colIdx)/sizeof(colIdx[0]));
    CSR_Handle *A = CSR_Create(m, m, rowPtr, colIdx, values);
    CSR_PrintInDense(A);
    printf("original banwdith: %d\n", CSR_GetBandwidth(A));

    printf("RCM permutation\n");
    int perm[m], inversePerm[m];
    CSR_GetRCMPemutation(A, perm, inversePerm);
    for (int i = 0; i < m; ++i) {
      printf("%d ", inversePerm[i]);
    }
    printf("\n");

    int rowPtr2[m + 1];
    int colIdx2[nnz];
    double values2[nnz];
    CSR_Handle *A2 = CSR_Create(m, m, rowPtr2, colIdx2, values2);

    CSR_Permute(A, A2, perm, inversePerm);
    CSR_PrintInDense(A2);
    printf("RCM permuted bandwidth: %d\n\n", CSR_GetBandwidth(A2));

    CSR_Destroy(A2);
    CSR_Destroy(A);
  }
#endif

  return 0;
}
