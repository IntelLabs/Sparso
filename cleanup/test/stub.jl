include("../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools.OptFramework

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_REORDER_WHEN_BENEFICIAL)

SparseAccelerator.dprintln(0, 1, "should not print \n as the \n level is too low")
SparseAccelerator.dprintln(1, 1, "1 tab")
SparseAccelerator.dprintln(1, 2, "2 tabs", "\nOK")
SparseAccelerator.dprintln(1, 3, "3 tabs", "\nOK")
flush(STDOUT)

function foo(A, x)
    A*x
end

m = 10
A = sprand(m, m, 0.1)
x = repmat([1/m], m)
@acc foo(A, x)
