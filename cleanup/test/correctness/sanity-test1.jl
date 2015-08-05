include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_VERBOSE)

SparseAccelerator.dprintln(0, 1, "should not print \n as the \n level is too low")
SparseAccelerator.dprintln(1, 0, "Testing tabbed display:")
SparseAccelerator.dprintln(1, 1, "1 tab")
SparseAccelerator.dprintln(1, 1, "OK")
SparseAccelerator.dprintln(1, 2, "2 tabs", "\nOK")
SparseAccelerator.dprintln(1, 3, "3 tabs", "\nOK")
flush(STDOUT)

