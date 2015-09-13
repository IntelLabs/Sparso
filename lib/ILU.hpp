#include "SpMP/CSR.hpp"
#include "SpMP/LevelSchedule.hpp"

namespace SpMP
{

// LU matrix that has the same sparsity structure as A is created.
// We only populate the value part of LU, and it shares rowptr and
// colidx with A.
// L can be constructed by tril(LU, -1) + speye(size(LU, 1)), so
// L is a unit lower triangular matrix
// U can be constructed by triu(LU)
void ilu0(double *lu, const CSR& A, const LevelSchedule& schedule);

// A version does the L and U splitting as well
void ilu0(CSR& L, CSR& U, const CSR& A, const LevelSchedule& schedule);

} // namespace SpMP
