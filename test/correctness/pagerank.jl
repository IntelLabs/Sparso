include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REORDER)

function pagerank(A, p, r) # p: initial rank, r: damping factor
  for i = 1:100
    p = (1-r) *A * p + r
  end
  p
end

A = matrix_market_read(ARGS[1], true, true)
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
r = 0.15

d = max(convert(Array{eltype(A),1}, vec(sum(A, 2))), 1) # num of neighbors
A = scale(A,1./d)

println("Original: ")
x = pagerank(A, p, r)
println("\tOriginal sum of x=", sum(x))

println("\nAccelerated: ")
@acc x= pagerank(A, p, r)
println("\tAccelerated sum of x=", sum(x))