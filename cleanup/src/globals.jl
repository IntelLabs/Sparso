# This file contains all the global constants, variables, routines.

typealias BasicBlock CompilerTools.CFGs.BasicBlock
typealias Statement  CompilerTools.CFGs.TopLevelStatement
typealias Liveness   CompilerTools.LivenessAnalysis.BlockLiveness
typealias CFG        CompilerTools.CFGs.CFG
typealias DomLoops   CompilerTools.Loops.DomLoops

# Options controlling debugging, performance (library choice, cost model), etc.
@doc """ Enable Sparse Accelerator """
const SA_ENABLE = 1

@doc """ Print verbose dump """
const SA_VERBOSE = 2

@doc """ Use Jula's default sparse matrix functions """
const SA_USE_JULIA = 8

@doc """ Use MKL sparse matrix functions """
const SA_USE_MKL = 16

@doc """ 
Use Sparse Matrix Pre-processing library (SpMP) functions. SPMP is a 
high-performance parallel implementation of BFS/RCM reordering, 
Gauss-Seidel smoother, sparse triangular solver, etc. 
"""
const SA_USE_SPMP = 32

@doc """ Enable reordering only when it is potentially beneficial """
const SA_REORDER_WHEN_BENEFICIAL = 64

# The internal booleans corresponding to the above options
sparse_acc_enabled      = false
verbosity               = 0
use_Julia               = false
use_MKL                 = false
use_SPMP                = true 
reorder_when_beneficial = false

@doc """ 
Set options for SparseAccelerator. The arguments can be any one or more 
of the following: SA_VERBOSE, SA_USE_JULIA, SA_USE_MKL, SA_USE_SPMP, 
SA_REORDER_WHEN_BENEFICIAL. They can appear in any order, except that 
SA_USE_JULIA, SA_USE_MKL and SA_USE_SPMP are exclusive with each other, and the
last one of them wins. 
"""
function set_options(args...)
    for arg in args
        if arg == SA_ENABLE
            if !sparse_acc_enabled
                # Insert sparse accelerator as 1 pass into the optimization framework
                sparse_accelerator_pass = OptFramework.optPass(SparseAccelerator.entry, true)
                OptFramework.setOptPasses([sparse_accelerator_pass])
            end
            global sparse_acc_enabled = true
        elseif arg == SA_VERBOSE 
            global verbosity = 1
            #OptFramework.set_debug_level(3)
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

@doc """
Build a dictionary for the symbols, including GenSym's, to store the type info 
of the typed func_ast. Note that all the symbols' scope is function-wise, even 
though in the source code, they might appear to have local or nested scopes: 
symbol renaming seems to have been done to make them function-wise.
"""
function build_symbol_dictionary(func_ast :: Expr)
    symbol_info = Dict{Union(Symbol, Integer), Type}()
    
    # Record Symbols' types
    assert(func_ast.head == :lambda)
    lambda = lambdaExprToLambdaInfo(func_ast)
    for i in lambda.var_defs
        symbol_info[i[2].name] = i[2].typ
    end
    
    # Record GenSym's types.
    for id = 1:length(lambda.gen_sym_typs)
        # Note that a GenSym id starts from 0
        symbol_info[id - 1] = lambda.gen_sym_typs[id] 
    end

    symbol_info
end

abstract Action

@doc """ Insert new statements to a basic block """
type InsertToBB <: Action
    new_stmts   :: Vector{Statement} 
    bb          :: BasicBlock
    at          :: Int 
end

@doc """ Insert new statements before a basic block """
type InsertBeforeBB <: Action
    new_stmts    :: Vector{Statement} 
    bb           :: BasicBlock
    outside_loop :: Bool # If bb is a loop header, outside_loop == true/false 
                         # would make the statements inserted outside/inside 
                         # the loop.
end

@doc """ Insert new statements on an edge """
type InsertOnEdge <: Action
    new_stmts    :: Vector{Statement} 
    predecessor  :: BasicBlock
    successor    :: BasicBlock
end

@doc """ 
All the analyses, including reordering, reusing, call replacement, etc.
"""
function analyses(
    func_ast    :: Expr, 
    symbol_info :: Dict{Union(Symbol,Integer), Type}, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops)
    
    actions = Vector{Action}()
    actions = reordering(actions, func_ast, symbol_info, liveness, cfg, loop_info)
    actions = context_info_discovery(actions, func_ast, symbol_info, liveness, 
                                     cfg, loop_info)
    actions
end

@doc """ 
Entry of SparseAccelerator. 
"""
function entry(func_ast :: Expr, func_arg_types :: Tuple, func_args)
    dprintln(1, 0, "******************************* SparseAccelerator ******************************")
    dprintln(1, 0, "Signature:")
    for i = 1 : length(func_args)
        if i == 1
            dprint(1, 1, "(", func_args[i], "::", func_arg_types[i]) 
        else
            dprint(1, 0, ", ", func_args[i], "::", func_arg_types[i])
        end 
    end
    dprintln(1, 0, ")\n\nAST:")
    dprintln(1, 1, func_ast)

    new_ast = nothing
    try
        # Build the common facilities: symbols' type dictionary, liveness
        # info, and control flow graph.
        symbol_info = build_symbol_dictionary(func_ast)
        liveness    = LivenessAnalysis.from_expr(func_ast)
        cfg         = liveness.cfg
        loop_info   = Loops.compute_dom_loops(cfg)

        # Do all analyses, and put their intended transformation code sequence
        # into a list. Then transform the code with the list.
        actions = analyses(func_ast, symbol_info, liveness, cfg, loop_info)
        new_ast = code_transformation(actions, func_ast, symbol_info, liveness, cfg, loop_info)

        dprintln(1, 0, "\nNew AST:")
        dprintln(1, 1, new_ast)
    catch ex
        dprintln(1, 0, "Exception! Sparse Accelerator skips optimizing the call.")
        dprintln(1, 1, ex)

        # Return the original AST without any change
        new_ast = func_ast
    finally
        dprintln(1, 0, "********************************************************************************")
        return new_ast
    end
end