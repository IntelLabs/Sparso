include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

function WAXPBY_test(x, alpha, p)
    x -= alpha * p
end

m = 10
x = repmat([1/m], m)
p = copy(x)
alpha = 0.1
@acc WAXPBY_test(x, alpha, p)
