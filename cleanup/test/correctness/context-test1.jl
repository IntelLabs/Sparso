include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)
CompilerTools.LivenessAnalysis.set_debug_level(6)

include("./pcg-symgs.jl")
include("utils.jl")

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1], true, true)
    m = size(A, 1)
end
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-10
maxiter = 1000

println("Original: ")
#x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
#println("\tsum of x=", sum(x))
#println("\tk=", k)
#println("\trel_err=", rel_err)

println("\nAccelerated: ")
x = zeros(Float64, m)
@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)