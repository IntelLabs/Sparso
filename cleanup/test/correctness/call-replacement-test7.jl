include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

function WAXPBY_test(p, beta, r)
    p = r + beta * p
end

m = 10
p = repmat([1/m], m)
r = copy(p)
beta = 0.1
@acc WAXPBY_test(p, beta, r)
