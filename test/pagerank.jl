#TODO: remove this include once the package is ready
include("../src/OptFramework.jl")
include("../src/SparseAccelerator.jl")
include("../src/sparse-analyze.jl")

using OptFramework
using MatrixMarket

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])

include("./mmread.jl")

function pagerank(A, p, r) # p: initial rank, r: damping factor
  tic()
  d = max(vec(sum(A, 2)), 1) # num of neighbors
  time1 = time()
  for i = 1:100
    q = p./d

#    tic()
    Aq = A*q
#    t += toq()

    p2 = r + (1-r)*Aq
#    println("$i: $(norm(p - p2)/norm(p))") # print out convergence

    p = p2
  end
  time2 = time()
  println("Time of original loop= ", time2 - time1, " seconds")  
  toc()
end

A = mmread(ASCIIString(ARGS[1]))
A = spones(A)

# unfortunately, colidx and rowptr in A are of Int64[], while our library assumes Cint[] (32 bit)
# convert to 32 bit first.
A1  = convert(SparseMatrixCSC{Cdouble, Cint}, A) 
A   = A1

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


