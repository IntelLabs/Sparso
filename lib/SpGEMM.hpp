#include "SpMP/CSR.hpp"

namespace SpMP
{

// C = A*diag(d)*B
CSR *inspectADB(const CSR *A, const CSR *B);
void adb(CSR *C, const CSR *A, const CSR *B, const double *d);

} // namespace SpMP
