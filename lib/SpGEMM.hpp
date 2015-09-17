#include "SpMP/CSR.hpp"

namespace SpMP
{

// C = A*diag(d)*B
void inspectADB(CSR *C, const CSR *A, const CSR *B);
CSR *inspectADB(const CSR *A, const CSR *B);
void adb(CSR *C, const CSR *A, const CSR *B, const double *d);
// C = A*B
CSR *SpGEMMWithEps(CSR *A, CSR *B, double eps);
// C = alpha*A + B
CSR *SpAdd(double alpha, CSR *A, double beta, CSR *B);

} // namespace SpMP
