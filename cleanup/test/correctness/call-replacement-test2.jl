include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

function SpMV_test(y, A, x)
    y = A * x
end

m = 10
A = sprand(m, m, 0.1)
x = repmat([1/m], m)
y = copy(x)
@acc SpMV_test(y, A, x)
