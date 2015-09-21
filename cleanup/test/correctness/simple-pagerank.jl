# This file is not complete. It is intended to be included into other file 
# to form a complete test

function pagerank(A, p, r) # p: initial rank, r: damping factor
  for i = 1 : 2 #100
    # SparseAccelerator.SpMV!(p, 1 - r, A, p, 0, p, r) # manual
    println("itr: ", i, " sumA=", sum(A), " sump=", sum(p),  " r=", r)
    p = (1-r) *A * p + r
    println("\tsump=", sum(p))
  end
  p
end