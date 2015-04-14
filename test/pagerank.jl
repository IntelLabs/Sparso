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

A = MatrixMarket.mmread(ARGS[1])
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
r = 0.15


tests = 5

println("**** original pagerank perf")
for i = 1 : tests
    p = repmat([1/m], m)
    pagerank(A, p, r)
end

println("**** accelerated pagerank perf")
for i = 1 : tests
    p = repmat([1/m], m)
    @acc pagerank(A, p, r)
end


