using CompilerTools
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using MatrixMarket
#TODO: remove this include once the package is ready
include("../src/SparseAccelerator.jl")

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])

function pagerank(A, p, r) # p: initial rank, r: damping factor
  repeat = 100
  Ap = copy(p)

  spmv_time = 0.0

  time1 = time()
  for i = 1:repeat
    spmv_time -= time()
    SparseAccelerator.SpMV!(Ap, 1 - r, A, p, 0, p, r) # manual
    spmv_time += time()

    if i == repeat
      println("Error = $(norm(p - Ap)/norm(p))") # print out convergence
    end

    temp = p
    p = Ap
    Ap = temp
  end
  time2 = time()

  original_loop_exec_time = time2 - time1
  println("Time of original loop= ", time2 - time1, " seconds")  
  println("SpMV time = $spmv_time, BW (GB/s) = $(nnz(A)*12.*100/spmv_time/1e9)\n")
  original_loop_exec_time
end

# This is for performance study only
function pagerank_reordered(A, p, r) # p: initial rank, r: damping factor
  spmv_time = 0
  d = max(vec(sum(A, 2)), 1) # num of neighbors

  time1 = time()
  
  __P_51227 = Array(Cint,size(A,2))
  __Pprime_51228 = Array(Cint,size(A,2))
  __A_51229 = SparseMatrixCSC(A.m,A.n,Array(Cint,size(A.colptr,1)),Array(Cint,size(A.rowval,1)),Array(Cdouble,size(A.nzval,1)))
  CSR_ReorderMatrix(A,__A_51229,__P_51227,__Pprime_51228,true,true,true)
  A = __A_51229

#  time2 = time()
  
  __p_51231 = Array(Cdouble,size(p,1))
  reorderVector(p,__p_51231,__P_51227)
  p = __p_51231
  
#  time3 = time()
  
  __d_51231 = Array(Cdouble,size(p,1))
  reorderVector(p,__d_51231,__P_51227)
  d = __d_51231
  
#  time4 = time()

#  println("**** Entering pagerank loop")
  for i = 1:100
    q = p./d
    time5 = time()
    Aq = A*q
    spmv_time += time() - time5
    p2 = r + (1-r)*Aq
    p = p2
  end
#  println("Exit pagerank loop *****")
  
  time5 = time()
  println("Time of original loop= ", time5 - time1, " seconds") 
  println("SpMV BW (GB/s) = ", nnz(A)*12.*100/spmv_time/1e9)
#  println("Breakdown:");
#  println("\tReorderMatrix A: ", time2 - time1, " seconds") 
#  println("\tReorderVector p: ", time3 - time2, " seconds") 
#  println("\tReorderVector d: ", time4 - time3, " seconds") 
#  println("\tLoop: ", time5 - time4, " seconds") 
end

A = MatrixMarket.mmread(ARGS[1])
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
r = 0.15

d = max(convert(Array{eltype(A),1}, vec(sum(A, 2))), 1) # num of neighbors
A = scale(A,1./d)

tests = 8
times = Float64[]

for lib in [SparseAccelerator.PCL_LIB]
    SparseAccelerator.use_lib(lib)
    
    original_time = 0.0
    for i = 1 : tests
        p = repmat([1/m], m)
        original_loop_exec_time = pagerank(A, p, r)
        original_time += original_loop_exec_time
    end
    
    accelerated_time = 0.0
    reorder_matrix_stats = Float64[0.0, 0.0, 0.0]
    ccall((:CSR_Statistics, "../lib/libcsr.so"), Void, (Ptr{Cdouble},), pointer(reorder_matrix_stats))
    for i = 1 : tests
        p = repmat([1/m], m)
        @acc original_loop_exec_time = pagerank(A, p, r)
    #    pagerank_reordered(A, p, r)
        accelerated_time += original_loop_exec_time
    end
        
    push!(times, original_time / tests)
    push!(times, accelerated_time / tests)
    for i = 1 : 3
        push!(times, reorder_matrix_stats[i]/tests)
    end
end

file = open("pagerank.timing", "a")
@printf(file, "%s", ARGS[1])
for i = 1 : size(times, 1)
    @printf(file, " %f", times[i])
end
@printf(file, "%s", "\n")
close(file)
