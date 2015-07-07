/* 
*   Matrix Market I/O library for ANSI C
*
*   See http://math.nist.gov/MatrixMarket for details.
*
*
*/

# ifndef MM_IO_H
# define MM_IO_H

#include "SpMP/mm_io.h"

#ifdef __cplusplus
extern "C" {
#endif

#define T double
// Julia should call the following two function in order.
// Between the two calls, the CSR array space must be allocated.
void load_matrix_market_step1 (char *file, int *sizes);
void load_matrix_market_step2 (char *file, T *a, int *j, int *i, int *sizes, bool one_based_CSR);

// C can directly call this once
void load_matrix_market (char *file, T **a, int **aj, int **ai, int *is_symmetric, int *am, int *an, int *annz, bool one_based_CSR = false);
#ifdef __cplusplus
}
#endif

# endif
