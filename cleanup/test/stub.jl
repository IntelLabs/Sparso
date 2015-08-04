include("../SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_VERBOSE, SA_USE_SPMP, SA_REORDER_WHEN_BENEFICIAL)

SparseAccelerator.dprint(0, 1, "should not print \n as the \n level is too low")
SparseAccelerator.dprint(1, 2, "should print \n and indent \n 8 blank spaces\n",
    "to be\n", "OK\n")

