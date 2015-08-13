# This file is not complete. It is intended to be included into other file 
# to form a complete test

include("utils.jl")

function pagerank(A, p, r) # p: initial rank, r: damping factor
  for i = 1:repeat
    # SparseAccelerator.SpMV!(p, 1 - r, A, p, 0, p, r) # manual
    p = (1-r) *A * p + r
  end
  p
end

m = 10
A = generate_symmetric_sparse_matrix(m)
p = repmat([1/m], m)
r = 0.15

