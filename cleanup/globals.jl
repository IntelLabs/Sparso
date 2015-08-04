# This file contains all the global constants, variables, routines.

typealias BasicBlock CompilerTools.CFGs.BasicBlock
typealias Statement  CompilerTools.CFGs.TopLevelStatement

# Options controlling debugging, performance (library choice, cost model), etc.
@doc """ Print verbose dump """
const SA_VERBOSE = 1

@doc """ Use Jula's default sparse matrix functions """
const SA_USE_JULIA = 8

@doc """ Use MKL sparse matrix functions """
const SA_USE_MKL = 16

@doc """ Use Sparse Matrix Pre-processing library (SpMP) functions. SPMP is a 
high-performance parallel implementation of BFS/RCM reordering, 
Gauss-Seidel smoother, sparse triangular solver, etc. """
const SA_USE_SPMP = 32

@doc """ Enable reordering only when it is potentially beneficial """
const SA_REORDER_WHEN_BENEFICIAL = 64

# The internal booleans corresponding to the above options
verbosity               = 0
use_Julia               = false
use_MKL                 = false
use_SPMP                = true 
reorder_when_beneficial = false

@doc """ Set options for SparseAccelerator. The arguments can be any one or more 
of the following: SA_VERBOSE, SA_USE_JULIA, SA_USE_MKL, SA_USE_SPMP, 
SA_REORDER_WHEN_BENEFICIAL. They can appear in any order, except that 
SA_USE_JULIA, SA_USE_MKL and SA_USE_SPMP are exclusive with each other, and the
last one of them wins. 
"""
function set_options(args...)
    for arg in args
        if arg == SA_VERBOSE 
            global verbosity = 1;
        elseif arg == SA_USE_JULIA 
            global use_Julia = true; global use_MKL = false; global use_SPMP = false
        elseif arg == SA_USE_MKL 
            global use_Julia = false; global use_MKL = true; global use_SPMP = false
        elseif arg == SA_USE_SPMP 
            global use_Julia = false; global use_MKL = false; global use_SPMP = true
        elseif arg == SA_REORDER_WHEN_BENEFICIAL 
            global reorder_when_beneficial = true
        else
            # TODO: print usage info
        end
    end
end

# Create a path to libcsr. This is a CSR(Compressed Sparse Row format)-based
# interface to the SPMP library.
const libcsr = joinpath(dirname(@__FILE__), "..", "lib", "libcsr.so")

@doc """ Entry of SparseAccelerator. """
function entry(func_ast :: Expr, func_arg_types :: Tuple, func_args)
end

# Insert sparse accelerator as 1 pass into the optimization framework
sparse_accelerator_pass = OptFramework.optPass(SparseAccelerator.entry, true)
OptFramework.setOptPasses([sparse_accelerator_pass])
