# This file contains all the global constants, variables, routines.

typealias BasicBlock CompilerTools.CFGs.BasicBlock
typealias Statement  CompilerTools.CFGs.TopLevelStatement

# Options controlling debugging, performance (library choice, cost model), etc.
@doc """ Print verbose dump """
const SA_VERBOSE = 1

@doc """ Use Jula's default sparse matrix functions """
const SA_USE_JULIA = 2

@doc """ Use MKL sparse matrix functions """
const SA_USE_MKL = 3

@doc """ Use Sparse Matrix Pre-processing library (SpMP) functions. SPMP is a 
high-performance parallel implementation of BFS/RCM reordering, 
Gauss-Seidel smoother, sparse triangular solver, etc. """
const SA_USE_SPMP = 4

@doc """ Enable reordering only when it is potentially beneficial """
const SA_REORDER_WHEN_BENEFICIAL = 5

# The internal booleans corresponding to the above options
verbose                 = false
use_Julia               = false
use_MKL                 = false
use_SPMP                = true 
reorder_when_beneficial = false

function set_options(args...)
    for arg in args
        if arg == 1 verbose = true
        elseif arg == 2 use_Julia = true; use_MKL = false; use_SPMP = false
        elseif arg == 3 use_Julia = false; use_MKL = true; use_SPMP = false
        elseif arg == 4 use_Julia = false; use_MKL = false; use_SPMP = true
        elseif arg == 5 reorder_when_beneficial = true
        else
            # TODO: print usage info
        end
    end
end

# Create a path to libcsr. This is a CSR(Compressed Sparse Row format)-based
# interface to the SPMP library.
const libcsr = joinpath(dirname(@__FILE__), "..", "lib", "libcsr.so")


# A debug print routine.
function dprint(level, msgs...)
    if(DEBUG_LVL >= level)
        print(msgs...)
    end
end

# A debug print routine.
function dprintln(level, msgs...)
    if(DEBUG_LVL >= level)
        println(msgs...)
    end
end


