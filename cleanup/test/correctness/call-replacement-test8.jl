include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

include("./simple-pagerank.jl")
include("utils.jl")

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1])
    m = size(A, 1)
end
p = repmat([1/m], m)
r = 0.15

println("Original: ")
p = pagerank(A, p, r)
println("\tsum of p=", sum(p))

println("\nAccelerated: ")
p = repmat([1/m], m)
@acc p = pagerank(A, p, r)
println("\tsum of p=", sum(p))

@acc pagerank(A, p, r)