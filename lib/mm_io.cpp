# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>

# include "mm_io.h"

/******************************************************************************/

int mm_is_valid ( MM_typecode matcode )

/******************************************************************************/
/*
  Purpose:

    MM_IS_VALID checks whether the MM header information is valid.

  Modified:

    31 October 2008

  Parameters:

    Input, MM_typecode MATCODE, the header information.

    Output, int MM_IS_VALID, is TRUE if the matrix code is valid.
*/
{
  if ( !mm_is_matrix ( matcode ) ) 
  {
    return 0;
  }

  if ( mm_is_dense ( matcode ) && mm_is_pattern ( matcode ) ) 
  {
    return 0;
  }

  if ( mm_is_real ( matcode ) && mm_is_hermitian ( matcode ) ) 
  {
    return 0;
  }

  if ( mm_is_pattern ( matcode ) && 
      ( mm_is_hermitian ( matcode ) || mm_is_skew ( matcode ) ) ) 
  {
    return 0;
  }
  return 1;
}
/******************************************************************************/

int mm_read_banner ( FILE *f, MM_typecode *matcode )

/******************************************************************************/
/*
  Purpose:

    MM_READ_BANNER reads the header line of an MM file.

  Modified:

    31 October 2008

  Parameters:

    Input, FILE *F, a pointer to the input file.

    Output, MM_typecode *MATCODE, the header information.
*/
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH]; 
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    mm_clear_typecode(matcode);  

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL) 
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, 
        storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p=mtx; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
    for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);
/* 
  check for banner 
*/
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

/* 
  first field should be "mtx" 
*/
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return  MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);

/* 
  second field describes whether this is a sparse matrix (in coordinate
  storgae) or a dense array 
*/
    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else
    if (strcmp(crd, MM_DENSE_STR) == 0)
            mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    

/*
  third field 
*/

    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else
    if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else
    if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    
/*
  fourth field 
*/

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else
    if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else
    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
        

    return 0;
}
/******************************************************************************/

int mm_read_mtx_array_size ( FILE *f, int *M, int *N )

/******************************************************************************/
/*
  Purpose:

    MM_READ_MTX_ARRAY_SIZE reads the size line of an MM array file.

  Modified:

    03 November 2008

  Parameters:

    Input, FILE *F, a pointer to the input file.

    Output, int *M, the number of rows, as read from the file.

    Output, int *N, the number of columns, as read from the file.

    Output, MM_READ_MTX_ARRAY_SIZE, an error flag.
    0, no error.
*/
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;
/* 
  set return null parameter values, in case we exit with errors 
*/
  *M = 0;
  *N = 0;
/* 
  now continue scanning until you reach the end-of-comments 
*/
  do 
  {
    if ( fgets ( line, MM_MAX_LINE_LENGTH, f ) == NULL ) 
    {
      return MM_PREMATURE_EOF;
    }
  } while ( line[0] == '%' );
/* 
  line[] is either blank or has M,N, nz 
*/
  if ( sscanf ( line, "%d %d", M, N ) == 2 )
  {
    return 0;
  }
  else
  {
    do
    { 
      num_items_read = fscanf ( f, "%d %d", M, N ); 
      if ( num_items_read == EOF ) 
      {
        return MM_PREMATURE_EOF;
      }
    }
    while ( num_items_read != 2 );
  }
  return 0;
}
/******************************************************************************/

int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I, int **J, 
        double **val, MM_typecode *matcode)

/******************************************************************************/
/*
  Purpose:

    MM_READ_MTX_CRD reads the values in an MM coordinate file.

  Discussion:

    This function allocates the storage for the arrays.

  Modified:

    31 October 2008

  Parameters:

*/
/*
    mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                        type code, e.g. 'MCRS'

                        if matrix is complex, values[] is of size 2*nz,
                            (nz pairs of real/imaginary values)
*/
{
    int ret_code;
    FILE *f;

    if (strcmp(fname, "stdin") == 0) f=stdin;
    else
    if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;


    if ((ret_code = mm_read_banner(f, matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) && 
            mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
        return ret_code;


    *I = (int *)  malloc(*nz * sizeof(int));
    *J = (int *)  malloc(*nz * sizeof(int));
    *val = NULL;

    if (mm_is_complex(*matcode))
    {
        *val = (double *) malloc(*nz * 2 * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }
    else if (mm_is_real(*matcode))
    {
        *val = (double *) malloc(*nz * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    else if (mm_is_pattern(*matcode))
    {
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    if (f != stdin) fclose(f);
    return 0;
}
/******************************************************************************/

int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
        double val[], MM_typecode matcode)

/******************************************************************************/
/*
  Purpose:

    MM_READ_MTX_CRD_DATA reads the values in an MM coordinate file.

  Discussion:

    This function assumes the array storage has already been allocated.

  Modified:

    31 October 2008

  Parameters:

    Input, FILE *F, a pointer to the input file.
*/
{
    int i;
    if (mm_is_complex(matcode))
    {
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
                != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
                != 3) return MM_PREMATURE_EOF;

        }
    }

    else if (mm_is_pattern(matcode))
    {
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d", &I[i], &J[i])
                != 2) return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}
/******************************************************************************/

int mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
        double *real, double *imag, MM_typecode matcode)

/******************************************************************************/
/*
  Purpose:

    MM_READ_MTX_CRD_ENTRY ???

  Modified:

    31 October 2008

  Parameters:

    Input, FILE *F, a pointer to the input file.
*/
{
    if (mm_is_complex(matcode))
    {
            if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
                != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
            if (fscanf(f, "%d %d %lg\n", I, J, real)
                != 3) return MM_PREMATURE_EOF;
    }
    else if (mm_is_integer(matcode))
    {
            long long integer;
            if (fscanf(f, "%d %d %lld\n", I, J, &integer)
                != 3) return MM_PREMATURE_EOF;
            *real = integer;
    }
    else if (mm_is_pattern(matcode))
    {
            if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}
/******************************************************************************/

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz )

/******************************************************************************/
/*
  Purpose:

    MM_READ_MTX_CRD_SIZE reads the size line of an MM coordinate file.

  Modified:

    03 November 2008

  Parameters:

    Input, FILE *F, a pointer to the input file.

    Output, int *M, the number of rows, as read from the file.

    Output, int *N, the number of columns, as read from the file.

    Output, int *NZ, the number of nonzero values, as read from the file.

    Output, MM_READ_MTX_CRF_SIZE, an error flag.
    0, no error.
*/
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;

/* 
  set return null parameter values, in case we exit with errors 
*/
  *M = 0;
  *N = 0;
  *nz = 0;
/* 
  now continue scanning until you reach the end-of-comments 
*/
  do 
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
    {
      return MM_PREMATURE_EOF;
    }
  } while (line[0] == '%');
/* 
  line[] is either blank or has M,N, nz 
*/
  if (sscanf(line, "%d %d %d", M, N, nz) == 3)
  {
    return 0;
  }
  else
  {
    do
    { 
      num_items_read = fscanf(f, "%d %d %d", M, N, nz); 
      if (num_items_read == EOF) 
      {
        return MM_PREMATURE_EOF;
      }
    } while (num_items_read != 3);
  }
  return 0;
}
/******************************************************************************/

int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
                double **val_, int **I_, int **J_)

/******************************************************************************/
/*
  Purpose:

    MM_READ_UNSYMMETRIC_SPARSE ???

  Modified:

    31 October 2008
*/
{
    FILE *f;
    MM_typecode matcode;
    int M, N, nz;
    int i;
    double *val;
    int *I, *J;
 
    if ((f = fopen(fname, "r")) == NULL)
            return -1;
 
 
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", fname);
        return -1;
    }
 
 
 
    if ( !(mm_is_real(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode)))
    {
        fprintf(stderr, "Sorry, this application does not support ");
        fprintf(stderr, "Market Market type: [%s]\n",
                mm_typecode_to_str(matcode));
        return -1;
    }
 
/* 
  find out size of sparse matrix: M, N, nz .... 
*/
 
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) !=0)
    {
        fprintf(stderr, "read_unsymmetric_sparse(): could not parse matrix size.\n");
        return -1;
    }
 
    *M_ = M;
    *N_ = N;
    *nz_ = nz;
/* 
  reseve memory for matrices 
*/ 
    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));
 
    *val_ = val;
    *I_ = I;
    *J_ = J;
/* 
  NOTE: when reading in doubles, ANSI C requires the use of the "l"
  specifier as in "%lg", "%lf", "%le", otherwise errors will occur
  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15) 
*/
 
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
/*
  Adjust from 1-based to 0-based.
*/
        I[i]--;
        J[i]--;
    }
    fclose(f);
 
    return 0;
}
/******************************************************************************/

char *mm_strdup ( const char *s )

/******************************************************************************/
/*
  Purpose:

    MM_STRDUP creates a new copy of a string.

  Discussion:

    This is a common routine, but not part of ANSI C, so it is included here.  
    It is used by mm_typecode_to_str().

  Modified:

    31 October 2008
*/
{
  int len = strlen ( s );
  char *s2 = (char *) malloc((len+1)*sizeof(char));
  return strcpy(s2, s);
}
/******************************************************************************/

char *mm_typecode_to_str ( MM_typecode matcode )

/******************************************************************************/
/*
  Purpose:

    MM_TYPECODE_TO_STR converts the internal typecode to an MM header string.

  Modified:

    31 October 2008
*/
{
    char buffer[MM_MAX_LINE_LENGTH];
    const char *types[4];
    char *mm_strdup(const char *);
    int error =0;

/* 
  check for MTX type 
*/
    if (mm_is_matrix(matcode)) 
        types[0] = MM_MTX_STR;
    else
        error=1;

/* 
  check for CRD or ARR matrix 
*/
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

/* 
  check for element data type 
*/
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;
/* 
  check for symmetry type 
*/
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else 
    if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else 
    if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return mm_strdup(buffer);

}
/******************************************************************************/

int mm_write_banner ( FILE *f, MM_typecode matcode )

/******************************************************************************/
/*
  Purpose:

    MM_WRITE_BANNER writes the header line of an MM file.

  Modified:

    31 October 2008

  Parameters:

    Input, FILE *F, a pointer to the output file.
*/
{
    char *str = mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    free(str);
    if (ret_code !=2 )
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}
/******************************************************************************/

int mm_write_mtx_array_size(FILE *f, int M, int N)

/******************************************************************************/
/*
  Purpose:

    MM_WRITE_MTX_ARRAY_SIZE writes the size line of an MM array file.

  Modified:

    31 October 2008
*/
{
    if (fprintf(f, "%d %d\n", M, N) != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}
/******************************************************************************/

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
        double val[], MM_typecode matcode)

/******************************************************************************/
/*
  Purpose:

    MM_WRITE_MTX_CRD writes an MM coordinate file.

  Modified:

    31 October 2008
*/
{
    FILE *f;
    int i;

    if (strcmp(fname, "stdout") == 0) 
        f = stdout;
    else
    if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;
/*
  print banner followed by typecode.
*/
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", mm_typecode_to_str(matcode));

/*
  print matrix sizes and nonzeros.
*/
    fprintf(f, "%d %d %d\n", M, N, nz);
/* 
  print values.
*/
    if (mm_is_pattern(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d\n", I[i], J[i]);
    else
    if (mm_is_real(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
    else
    if (mm_is_complex(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i], 
                        val[2*i+1]);
    else
    {
        if (f != stdout) fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

  if ( f !=stdout ) 
  {
    fclose(f);
  }
  return 0;
}
/******************************************************************************/

int mm_write_mtx_crd_size ( FILE *f, int M, int N, int nz )

/******************************************************************************/
/*
  Purpose:

    MM_WRITE_MTX_CRD_SIZE writes the size line of an MM coordinate file. 

  Modified:

    31 October 2008

  Parameters:

    Input, FILE *F, a pointer to the output file.
*/
{
  if ( fprintf ( f, "%d %d %d\n", M, N, nz ) != 3 )
  {
    return MM_COULD_NOT_WRITE_FILE;
  }
  else 
  {
    return 0;
  }
}

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

#define RAND01() randfp(0.0, 1.0)

void print_some_COO_values(int nnz, double *values, int *rowidx, int *colidx)
{
  fflush(stdout);
  int distance = nnz/100; // print about 100 elements
  printf("COO values: \n");
  for (int i = 0; i < nnz; i++) {
    if ((i % distance) == 0 || i < 10 || nnz - i < 10) {
      printf("%d %d %.11f\n", rowidx[i], colidx[i], values[i]);
    }
  }
  fflush(stdout);
}

// COO: true: build a COO array, stored statically inside. 
//      false: convert the COO array to the CSR array (a, aj, ai). Then free the COO array space
// This function should be called twice, first with COO as true, next as false.
// Between the two calls, the CSR array space must be allocated.
void load_matrix_market_step (char *file, T *a, int *j, int *i, int *sizes, bool COO, bool one_based_CSR)
{
    // The COO array
    static int *colidx;
    static int *rowidx;
    static T *values;
        
    if (COO) {
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
                values[count] = values[count -1]; // it must be symmetric no matter it is pattern or not
                count++;
                trc++;
            }
    
        }
        assert (lines <= nnz);
        nnz = count;

        //print_some_COO_values(nnz, values, rowidx, colidx);
    
        sizes[0] = (int)symm;
        sizes[1] = m;
        sizes[2] = n;
        sizes[3] = nnz;
        
        assert(!zero_based_ind);
        const char *tag=(symm?"symmetric":"general");
        printf("%s:::%s m=%d n=%d nnz=%d\n", file, tag, m, n, nnz);
        fclose(fp);
    } else {
        // Convert COO to CSR, then free COO
        int m = sizes[1];
        int nnz = sizes[3];

        printf("************* step 2, some COO values: m =%d, nnz=%d\n", m, nnz);        
        //print_some_COO_values(nnz, values, rowidx, colidx);
        
        coo2csr(m, nnz, values, rowidx, colidx, a, j, i);
        if (one_based_CSR) {
            set_1based_ind(i, j, m, nnz);
        }
        
        free (colidx);
        free (rowidx);
        free (values);
    }
}

// The following two functions are for Julia
void load_matrix_market_step1 (char *file, int *sizes)
{
    load_matrix_market_step(file, NULL, NULL, NULL, sizes, true, true /*useless here*/);
}

void load_matrix_market_step2 (char *file, T *a, int *j, int *i, int *sizes, bool one_based_CSR)
{
    load_matrix_market_step(file, a, j, i, sizes, false, one_based_CSR);
}

void load_matrix_market (char *file, T **a, int **aj, int **ai, int *is_symmetric, int *am, int *an, int *annz, bool one_based_CSR)
{
    // Read into a COO array (stored statically insded mm_read)
    int sizes[] = {0, 0, 0, 0};
    load_matrix_market_step1(file, sizes);
    
    *is_symmetric = sizes[0];
    *am = sizes[1];
    *an = sizes[2];
    *annz = sizes[3];
    
    // Allocate space for the CSR array
    *a = (T *)_mm_malloc(sizeof(T)*(*annz), 64);
    *aj = (int *)_mm_malloc(sizeof(int)*(*annz), 64);
    *ai = (int *)_mm_malloc(sizeof(int)*(*am+1), 64);

    // Convert COO to CSR
    load_matrix_market_step2(file, *a, *aj, *ai, sizes, one_based_CSR);
}

/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
