include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REORDER_WHEN_BENEFICIAL)

function foo(A, x)
    A*x
end

m = 10
A = sprand(m, m, 0.1)
x = repmat([1/m], m)
@acc foo(A, x)
