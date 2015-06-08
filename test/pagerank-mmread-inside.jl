using CompilerTools
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using MatrixMarket
#TODO: remove this include once the package is ready
include("../src/SparseAccelerator.jl")

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])
SparseAccelerator.set_debug_level(2)

function pagerank(ARGS)
  A = MatrixMarket.mmread(ARGS[1])
  A = spones(A)

  m = size(A, 1)
  
  # p: initial rank, r: damping factor
  p = repmat([1/m], m)
  r = 0.15

  d = max(convert(Array{eltype(A),1}, vec(sum(A, 2))), 1) # num of neighbors
  A = scale(A, 1./d)

  time1 = time()
  for i = 1:100
    q = p./d
    Aq = A*q
    p = r + (1-r)*Aq
  end
  time2 = time()

  original_loop_exec_time = time2 - time1
  println("Time of original loop= ", time2 - time1, " seconds")  
  original_loop_exec_time
end

tests = 1
times = Float64[]

for lib in [SparseAccelerator.PCL_LIB]
    SparseAccelerator.use_lib(lib)
    
    original_time = 0.0
    for i = 1 : tests
        original_loop_exec_time = 0 #pagerank(ARGS)
        original_time += original_loop_exec_time
    end
    
    accelerated_time = 0.0
    reorder_matrix_stats = Float64[0.0, 0.0, 0.0]
    ccall((:CSR_Statistics, "../lib/libcsr.so"), Void, (Ptr{Cdouble},), pointer(reorder_matrix_stats))
    for i = 1 : tests
        @acc original_loop_exec_time = pagerank(ARGS)
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
