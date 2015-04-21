#TODO: remove this include once the package is ready
include("../src/OptFramework.jl")
include("../src/SparseAccelerator.jl")
include("../src/sparse-analyze.jl")

using OptFramework
using MatrixMarket

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])

function pagerank(A, p, r) # p: initial rank, r: damping factor
  tic()
  # The following convert is needed so Julia doesn't give the result of "vec" to be of type Array{T,N}.
  # We are smarter here and convert to exactly the right type.  Without this convert, d and q will be
  # of a union type and SpMV using q won't be recognized as distributive.
  d = max(convert(Array{eltype(A),1}, vec(sum(A, 2))), 1) # num of neighbors
  time1 = time()
  for i = 1:100
    q = p./d

#    tic()
#    Aq = A*q
     Aq = SpMV(A, q)
#    t += toq()

    p2 = r + (1-r)*Aq
#    println("$i: $(norm(p - p2)/norm(p))") # print out convergence

    p = p2
  end
  time2 = time()
  println("Time of original loop= ", time2 - time1, " seconds")  
  toc()
end

# This is for performance study only
function pagerank_reordered(A, p, r) # p: initial rank, r: damping factor
  d = max(vec(sum(A, 2)), 1) # num of neighbors

  time1 = time()
  
  @time __P_51227 = Array(Cint,size(A,2))
  @time __Pprime_51228 = Array(Cint,size(A,2))
  @time __A_51229 = SparseMatrixCSC(A.m,A.n,Array(Cint,size(A.colptr,1)),Array(Cint,size(A.rowval,1)),Array(Cdouble,size(A.nzval,1)))
  @time CSR_ReorderMatrix(A,__A_51229,__P_51227,__Pprime_51228,true,true,true)
  @time A = __A_51229

  time2 = time()
  
  __p_51231 = Array(Cdouble,size(p,1))
  reorderVector(p,__p_51231,__P_51227)
  p = __p_51231
  
  time3 = time()
  
  __d_51231 = Array(Cdouble,size(p,1))
  reorderVector(p,__d_51231,__P_51227)
  d = __d_51231
  
  time4 = time()

  println("**** Entering pagerank loop")
  for i = 1:100
    q = p./d
    @time Aq = A*q
    p2 = r + (1-r)*Aq
    p = p2
  end
  println("Exit pagerank loop *****")
  
  time5 = time()
  println("Time of original loop= ", time5 - time1, " seconds") 
  println("Breakdown:");
  println("\tReorderMatrix A: ", time2 - time1, " seconds") 
  println("\tReorderVector p: ", time3 - time2, " seconds") 
  println("\tReorderVector d: ", time4 - time3, " seconds") 
  println("\tLoop: ", time5 - time4, " seconds") 
end

A = MatrixMarket.mmread(ARGS[1])
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
r = 0.15


tests = 1

println("**** original pagerank perf")
for i = 1 : tests
    p = repmat([1/m], m)
    pagerank(A, p, r)
end

println("**** accelerated pagerank perf")
for i = 1 : tests
    p = repmat([1/m], m)
#    @acc pagerank(A, p, r)
    pagerank_reordered(A, p, r)
end


