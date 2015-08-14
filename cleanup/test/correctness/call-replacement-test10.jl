include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

include("./pcg-jacobi.jl")
@acc pcg_jacobi(x, A, b, tol, maxiter)

