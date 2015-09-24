include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP)
include("utils.jl")

function pagerank(A, p, r) # p: initial rank, r: damping factor
  for i = 1 : 2 #100
    p = (1-r) *A * p + r
  end
  p
end

function pagerank_manual(A, p, r) # p: initial rank, r: damping factor
  for i = 1 : 2 #100
    SparseAccelerator.SpMV!(p, 1 - r, A, p, 0, p, r) # manual
  end
  p
end

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1], true, true)
    m = size(A, 1)
end
r = 0.15

println("Original: ")
p = repmat([1/m], m)
p = pagerank(A, p, r)
println("\tOriginal sum of p=", sum(p))

println("Manual: ")
p = repmat([1/m], m)
p = pagerank_manual(A, p, r)
println("\tManual sum of p=", sum(p))