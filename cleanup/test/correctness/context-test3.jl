include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using MatrixMarket

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)
CompilerTools.LivenessAnalysis.set_debug_level(6)

include("./ipm-ref.jl")
include("utils.jl")

if length(ARGS) == 0
    # min 2*x1 + x2 subject to x1 + x2 = 1, x1 >= 0, x2 >= 0
    # expected solution: x1 = 0, x2 = 1, obj = 1
    A = sparse([1 1])
    b = [ 1 ]'
    p = [ 2 1 ]'
else
    A = matrix_market_read(string(ARGS[1], "-A.mtx"))'
    b = vec(MatrixMarket.mmread(string(ARGS[1], "-b.mtx")))
    p = vec(MatrixMarket.mmread(string(ARGS[1], "-p.mtx")))
end

m = size(A, 1)
n = size(A, 2)
println("Problem size = [$m $n]")
#println("\tsum of A=", sum(A))
#println("\tsum of b=", sum(b))
#println("\tsum of p=", sum(p))

#println("Original: ")
#x = ipm_ref1(A, b, p)
#println("\tsum of x=", sum(x))

println("\nAccelerated: ")
@acc x = ipm_ref(A, b, p)
println("\tsum of x=", sum(x))






