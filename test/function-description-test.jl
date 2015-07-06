using CompilerTools
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using MatrixMarket

include("../Src/SparseAccelerator.jl")
include("./MatrixMarket.jl")
using SparseAccelerator

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])
SparseAccelerator.set_debug_level(2)

function test_SpMV(A, x)
    y = copy(x)
    for i = 1:10
        y += SparseAccelerator.SpMV(A, x)
    end
    return y
end

m = 10
A = sprand(m, m, 0.1)
# Make A symmetric, which is the assumption of reordering for now
for i = 1:m for j = 1:m A[i, j] = A[j, i] end end
x = repmat([1/m], m)
@acc test_SpMV(A, x)
