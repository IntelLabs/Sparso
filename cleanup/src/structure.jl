# Some patterns are for capturing structrures of matrices.

@doc """
Describe part (lower, upper, diagonal) of a matrix, of which the structure matters.
TODO: add symmetric properpty.
"""
type StructureProxy
    lower    :: Bool
    upper    :: Bool
    diagonal :: Bool
    proxy    :: Any
end

structure_proxies = Dict{Any, StructureProxy}()

@doc """" Set for a matrix of AST a structure proxy in the dictionary """
function set_structure_proxy(
    A         :: Any,
    structure :: StructureProxy
)
    if typeof(A) == SymbolNode
        structure_proxies[A.name] = structure
    else
        structure_proxies[A] = structure
    end
end

@doc """" Search for a matrix of AST a structure proxy in the dictionary """
function get_structure_proxy(
    A :: Any
)
    if typeof(A) == SymbolNode
        A = A.name
    end
    if !haskey(structure_proxies, A)
        return nothing
    end
    return structure_proxies[A]
end

@doc """ Post-processing function. Propagate lower part of a matrix. """
function propagate_lower_structure(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    proxy = ast.args[2]
    set_structure_proxy(ast, StructureProxy(true, false, false, proxy))
    return true
end

@doc """ Post-processing function. Propagate upper part of a matrix. """
function propagate_upper_structure(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    proxy = ast.args[2]
    set_structure_proxy(ast, StructureProxy(false, true, false, proxy))
    return true
end

@doc """ Post-processing function. Memoize a diagonal matrix. """
function memoize_diagonal_structure(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    set_structure_proxy(ast, StructureProxy(false, false, true, nothing))
    return true
end

@doc """ Check if arg2 is a sparse diagonal matrix """
function CS_spdiagm_times_any_check(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    A = ast.args[2]
    structure = get_structure_proxy(A) 
    if structure == nothing || !(structure.diagonal)
        return false
    end
    return true
end

@doc """ Post-processing function. Propagate the structure of the last arg. """
function propagate_last_structure(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    A = last(ast.args)
    structure = get_structure_proxy(A) 
    if structure != nothing
        if ast.head == :(=)
            set_structure_proxy(ast.args[1], structure_proxies[A])
        else
            set_structure_proxy(ast, structure_proxies[A])
        end
    end
    return true
end

@doc """ A function that will be called when no pattern could handle an AST. """
function last_resort(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    # TODO: put whatever you want as the last resort of the analysis here.
    return true
end

const CS_tril_pattern = ExprPattern(
    "CS_tril_pattern",
    (:call, GlobalRef(Main, :tril), SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    propagate_lower_structure,
    "",
    ""
)

const CS_triu_pattern = ExprPattern(
    "CS_triu_pattern",
    (:call, GlobalRef(Main, :triu), SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    propagate_upper_structure,
    "",
    ""
)

const CS_spdiagm_pattern = ExprPattern(
    "CS_spdiagm_pattern",
    (:call, GlobalRef(Main, :spdiagm), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    memoize_diagonal_structure,
    "",
    ""
)

const CS_spdiagm_times_any_pattern = ExprPattern(
    "CS_spdiagm_times_any_pattern",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    CS_spdiagm_times_any_check,
    (:NO_CHANGE, ),
    propagate_last_structure,
    "",
    ""
)

const CS_assign_pattern = ExprPattern(
    "CS_assign_pattern",
    (:(=), Any, Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    propagate_last_structure,
    "",
    ""
)

# This is the only pattern that will always be matched, justed based on its name.
# It should always be put as the last pattern, and when it is matched, 
# the last_resort() function will be called.
const CS_last_resort_pattern = ExprPattern(
    "CS_last_resort_pattern", # Useless
    (),                       # Useless
    (:NO_SUB_PATTERNS,),      # Useless
    last_resort,
    (:NO_CHANGE, ),           # Useless
    do_nothing,               # Useless
    "",                       # Useless
    ""                        # Useless
)

@doc """
Patterns used for discovering matrix structures. 
"""
CS_structure_propagation_patterns = [
    CS_tril_pattern,
    CS_triu_pattern,
    CS_spdiagm_pattern,
    CS_spdiagm_times_any_pattern,
    CS_assign_pattern,
    CS_last_resort_pattern
]

@doc """ 
Discover matrix structures in the whole function.
"""
function structure_discovery(
    symbol_info :: Sym2TypeMap, 
    cfg         :: CFG
)

    call_sites  = CallSites(Set{CallSite}(), WholeFunction(), symbol_info, Set{Sym}(),
                            CS_structure_propagation_patterns, Vector{Action}())
    for (bb_idx, bb) in cfg.basic_blocks
        for stmt in bb.statements
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end
            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
        end
    end
    
    dprintln(1, 0, "\nMatrix structures discovered:")
    dprintln(1, 1, "", structure_proxies)
end
