# This file contains all the global constants, variables, routines.

typealias BasicBlock      CompilerTools.CFGs.BasicBlock
typealias Statement       CompilerTools.CFGs.TopLevelStatement
typealias Liveness        CompilerTools.LivenessAnalysis.BlockLiveness
typealias CFG             CompilerTools.CFGs.CFG
typealias DomLoops        CompilerTools.Loops.DomLoops
typealias Loop            CompilerTools.Loops.Loop
typealias GenSymId        Int
typealias BasicBlockIndex Int
typealias StatementIndex  Int
typealias Sym             Union(Symbol, GenSymId) # A Symbol or a GenSym id.
typealias Sym2TypeMap     Dict{Sym, Type}

# Options controlling debugging, performance (library choice, cost model), etc.
@doc """ Enable Sparse Accelerator """
const SA_ENABLE = 0

@doc """ Print verbose dump """
const SA_VERBOSE = 8

@doc """ Use Jula's default sparse matrix functions """
const SA_USE_JULIA = 16

@doc """ Use MKL sparse matrix functions """
const SA_USE_MKL = 24

@doc """ 
Use Sparse Matrix Pre-processing library (SpMP) functions. SPMP is a 
high-performance parallel implementation of BFS/RCM reordering, 
Gauss-Seidel smoother, sparse triangular solver, etc. 
"""
const SA_USE_SPMP = 32

@doc """" 
Pattern match and replace the code that is functionally equivalent to SpMV, dot,
WAXPBY, etc. with calls to the corresponding SPMP library functions.
"""
const SA_REPLACE_CALLS = 40

@doc """ Enable reordering of arrays. """
const SA_REORDER = 48

@doc """ Enable reordering only when it is potentially beneficial """
const SA_REORDER_WHEN_BENEFICIAL = 56

@doc """ Enable context-sensitive optimization. """
const SA_CONTEXT = 64


# The internal booleans corresponding to the above options, and their default values
sparse_acc_enabled            = false
verbosity                     = 0
use_Julia                     = false
use_MKL                       = false
use_SPMP                      = true
replace_calls_enabled         = false
reorder_enabled               = false
reorder_when_beneficial       = false
context_sensitive_opt_enabled = false

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
        elseif arg == SA_REPLACE_CALLS
            global replace_calls_enabled = true
        elseif arg == SA_REORDER
            global reorder_enabled = true
        elseif arg == SA_REORDER_WHEN_BENEFICIAL 
            global reorder_when_beneficial = true
        elseif arg == SA_CONTEXT
            global context_sensitive_opt_enabled = true
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
    symbol_info = Sym2TypeMap()
    
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

@doc """ 
Determine the type of an AST node. A Symbol or GenSym gets a type from the 
symbol_info. An expression or SymbolNode gets the type stored by Julia type 
inference. All the other kinds of AST nodes resort to the default typeof(). 
"""
function type_of_ast_node(node, symbol_info :: Sym2TypeMap)
    local typ = typeof(node)
    if typ == Symbol
        # Use get() instead [] in case the key (like Symbol "*") does not exist
        # Return Nothing if no info found
        return get(symbol_info, node, Nothing)
    elseif typ == GenSym
        return get(symbol_info, node.id, Nothing)
    elseif typ == Expr || typ == SymbolNode
        return node.typ
    else
        return typ
    end
end

@doc """ Is the type a scalar (number), or an array? """
function is_number_or_array(typ :: Type)
    is_number = (typ <: Number)

    # We assume the user program does not use Range for any array
    # computation, although Range is a subtype of AbstractArray
    is_array  = (typ <: AbstractArray && !(typ <: Range))

    is_number, is_array
end

@doc """ Are the types scalars (numbers), or are some of them arrays? """
function are_numbers_or_arrays(result_type :: Type, arg_types :: Tuple{Type})
    all_numbers, some_arrays = is_number_or_array(result_type)
    for t in arg_types
        is_number, is_array = is_number_or_array(t)
        all_numbers         = all_numbers && is_number
        some_arrays         = some_arrays || is_array
    end
    all_numbers, some_arrays
end

@doc """ A module (or function)'s name string """
function module_or_function_name(arg)
    if typeof(arg) == Symbol
        return string(arg)
    elseif typeof(arg) == QuoteNode && typeof(arg.value) == Symbol
        return string(arg.value)
    else
        throw(UnhandledModuleOrFunctionName(arg))
    end
end

@doc """ 
From the given call arguments, figure out the module and function name 
"""
function resolve_call_names(call_args :: Vector)
    assert(length(call_args) > 0)
    module_name, function_name = "", ""
    if isdefined(:GlobalRef) && typeof(call_args[1]) == GlobalRef
        return string(call_args[1].mod), string(call_args[1].name)
    elseif typeof(call_args[1]) == Symbol 
        # Example: :*
        function_name = string(call_args[1])
    elseif isa(call_args[1], TopNode) && length(call_args) == 3
        # Example: top(getfield), SparseAccelerator,:SpMV
        if call_args[1] == TopNode(:getfield)
            module_name   = module_or_function_name(call_args[2])
            function_name = module_or_function_name(call_args[3])
        else
            function_name = module_or_function_name(call_args[1].name)
        end
    elseif  isa(call_args[1], Expr) &&
        call_args[1].head == :call
        # Example: (:call, top(getfield), SparseAccelerator,:SpMV)
        return resolve_call_names(call_args[1].args)
    end
    module_name, function_name
end

@doc """
Search in a database for information of the function specified by the module 
name, function name, and argument types.
"""
function look_for_function(
    database       :: Vector, 
    module_name    :: String, 
    function_name  :: String, 
    argument_types :: Tuple{Type}
)
    for item in database
        if module_name    == item.module_name   && 
           function_name  == item.function_name &&
           length(argument_types) == length(item.argument_types)
            found = true
            for i in 1:length(argument_types) 
               if !(argument_types[i] <: item.argument_types[i])
                   found = false
                   break
                end
            end
            if found
               return item
            end
        end
    end
    return nothing
end

@doc """
A call sites of interesting functions (like SpMV, triangular solver, etc.). From
the AST of the call, we may figure out the matrices in its arguments so that we 
may create a matrix knob for it. We may also create a function knob for the call
site, and may delete the knob later.
"""
type CallSite
    ast           :: Expr 
    matrices      :: Set{Sym} 
    fknob_creator :: Tuple{Symbol, String} # A library call to create a function knob for this call site
    fknob_deletor :: Tuple{Symbol, String} # A library call to delete the function knob for this call site
end

@doc """
Call sites of interesting functions (like SpMV, triangular solver, etc.). The
function's result and argument types are figured out with the help of symbol_info. 
"""
type CallSites
    sites       :: Set{CallSite}
    symbol_info :: Sym2TypeMap
end

abstract Action

@doc """ Insert new statements to a basic block """
type InsertBeforeStatement <: Action
    new_stmts   :: Vector{Statement} 
    bb          :: BasicBlock
    stmt        :: Statement 
end

@doc """ Insert new statements before a loop's head block """
type InsertBeforeLoopHead <: Action
    new_stmts    :: Vector{Statement} 
    loop         :: Loop # The loop.
    outside_loop :: Bool # Outside_loop == true/false would make the statements
                         # inserted outside/inside the loop.
end

@doc """ Insert new statements on an edge """
type InsertOnEdge <: Action
    new_stmts    :: Vector{Statement} 
    from_bb      :: BasicBlock
    to_bb        :: BasicBlock
end

@doc """ 
All the analyses, including reordering, reusing, call replacement, etc.
"""
function analyses(
    func_ast    :: Expr, 
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops)

    dprintln(1, 0, "\nAnalyzing ...", "\nCFG:")
    dprintln(1, 1, cfg)
        
    regions = region_formation(cfg, loop_info)
    actions = Vector{Action}()
    if reorder_enabled
        actions = reordering(actions, regions, func_ast, symbol_info, liveness, cfg, loop_info)
    end
    if context_sensitive_opt_enabled
        actions = context_info_discovery(actions, regions, func_ast, symbol_info, liveness, cfg)
    end
    actions
end

@doc """ 
Entry of SparseAccelerator. 
"""
function entry(func_ast :: Expr, func_arg_types :: Tuple, func_args)
    old_ast = copy(func_ast)
    new_ast = nothing
    try
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
    
        assert(func_ast.head == :lambda)
        
        # Build the common facilities: symbols' type dictionary, liveness
        # info, and control flow graph.
        symbol_info = build_symbol_dictionary(func_ast)
        liveness    = LivenessAnalysis.from_expr(func_ast)
        cfg         = liveness.cfg
        loop_info   = Loops.compute_dom_loops(cfg)

        dprintln(1, 0, "\nFunction body showing structures:")
        dsprintln(1, 1, symbol_info, func_ast)

        new_ast = func_ast
        
        if replace_calls_enabled
            replace_calls(symbol_info, cfg)
        end
        
        if reorder_enabled || context_sensitive_opt_enabled
            # Reordering and context-sensitive optimization: Do all analyses, and 
            # put their intended transformation code sequence into a list. Then 
            # transform the code with the list.
            actions = analyses(func_ast, symbol_info, liveness, cfg, loop_info)
            new_ast = code_transformation(actions, func_ast, cfg)
        end
        
        dprintln(1, 0, "\nNew AST:")
        dprintln(1, 1, new_ast)
    catch ex
        dprintln(1, 0, "Exception! Sparse Accelerator skips optimizing the call.")
        dprintln(1, 1, ex)

        # Return the original AST without any change
        new_ast = old_ast
    finally
        dprintln(1, 0, "********************************************************************************")
        return new_ast
    end
end