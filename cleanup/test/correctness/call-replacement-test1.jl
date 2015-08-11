include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REPLACE_CALLS)

function dot_test(x, y)
    dot(x, y)
end

m = 10
x = repmat([1/m], m)
y = copy(x)
@acc dot_test(x, y)
