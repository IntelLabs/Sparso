include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

function SpMV!_test(y, A, x)
    A_mul_B!(y, A, x)
end

m = 10
A = sprand(m, m, 0.1)
x = repmat([1/m], m)
y = copy(x)
@acc SpMV!_test(y, A, x)
