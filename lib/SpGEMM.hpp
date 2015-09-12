#include "SpMP/CSR.hpp"

namespace SpMP
{

// C = A*diag(d)*B
void inspectADB(CSR *C, const CSR *A, const CSR *B);
CSR *inspectADB(const CSR *A, const CSR *B);
void adb(CSR *C, const CSR *A, const CSR *B, const double *d);

} // namespace SpMP
