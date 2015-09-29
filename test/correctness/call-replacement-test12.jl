include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

include("./pcg.jl")

A       = matrix_market_read(ARGS[1], true, true)
m       = size(A, 1)
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-7
maxiter = 20000

M = A # Perfect


@acc pcg(x, A, b, M, tol, maxiter)

