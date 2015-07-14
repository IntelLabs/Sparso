include("../Src/SparseAccelerator.jl")
include("./MatrixMarket.jl")
using SparseAccelerator

using CompilerTools
using SparseAccelerator
using CompilerTools.OptFramework
using MatrixMarket
using Base.Test

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])
#CompilerTools.LivenessAnalysis.set_debug_level(3)
#SparseAccelerator.set_debug_level(3)
#OptFramework.set_debug_level(3)
SparseAccelerator.set_debug_level(2)

function pcg(x, A, b, M, tol, maxiter)
    tic()
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = M \ r
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p
#        Ap = SparseAccelerator.SpMV(A, p) # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end
        z = M \ r  
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    time2 = time()
    original_loop_exec_time = time2 - time1
    println("Time of original loop= ", time2 - time1, " seconds")
    toc()
    return x, k, rel_err, original_loop_exec_time
end

#matrices = [
#"hpcg_4",
#"hpcg_32",
#]

tol = 1e-10
maxiter = 1000

SparseAccelerator.use_lib(SparseAccelerator.PCL_LIB)

#for matrix in matrices
#  println(matrix)
#  A = MatrixMarket.mmread("data/$matrix.mtx")

#  A = MatrixMarket.mmread(ARGS[1])
A = sprand(10, 10, 1.0)
#  A = spones(A)

  N   = size(A, 1)
for i = 1:N for j=1:N A[i,j] = A[j,i] end end 
  
  b   = ones(Float64, N)

  tests = 1
  times = Float64[]

  println("\n\n**** original pcg perf")
  loop_exec_time = 0.0
  for i = 1:tests  
      M = A # Perfect
      x   = zeros(Float64, N)
      x, k, err, original_loop_exec_time = pcg(x, A, b, M, tol, maxiter)
      loop_exec_time += original_loop_exec_time
      println("Perfect preconditioner: $k iterations $err error $original_loop_exec_time sec.")
  end
  push!(times, loop_exec_time / tests)

  println("\n\n**** accelerated pcg perf")
  loop_exec_time = 0.0
  for i = 1:tests  
      M = A # Perfect
      x   = zeros(Float64, N)
      @acc result = pcg(x, A, b, M, tol, maxiter)
      x, k, err, original_loop_exec_time = result
      loop_exec_time += original_loop_exec_time
      println("Perfect preconditioner: $k iterations $err error $original_loop_exec_time sec.")
  end
  push!(times, loop_exec_time / tests)
  
  println()

file = open("pcg.timing", "a")
@printf(file, "%s", ARGS[1])
for i = 1 : size(times, 1)
    @printf(file, " %f",  times[i])
end
@printf(file, "%s", "\n")
close(file)