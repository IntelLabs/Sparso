# This file is not complete. It is intended to be included into other file 
# to form a complete test

function pagerank(A, p, r) # p: initial rank, r: damping factor
  for i = 1 : 100
    # SparseAccelerator.SpMV!(p, 1 - r, A, p, 0, p, r) # manual
    p = (1-r) *A * p + r
  end
  p
end