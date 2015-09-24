include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("./pcg-symgs.jl")
include("utils.jl")

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1], true, true)
    m = size(A, 1)
end
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-10
maxiter = 1000

println("Original: ")
x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tOriginal sum of x=", sum(x))
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

println("\nAccelerated: ")
x = zeros(Float64, m)
@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tAccelerated sum of x=", sum(x))
println("\tAccelerated k=", k)
println("\tAccelerated rel_err=", rel_err)