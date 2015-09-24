include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

include("./pcg.jl")
@acc pcg(x, A, b, M, tol, maxiter)

