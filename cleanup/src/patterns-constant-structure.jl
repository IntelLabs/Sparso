@doc """ Post-processing function. Propagate lower part of a matrix. """
function propagate_lower_structure(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple
)
    proxy = ast.args[2]
    set_structure_proxy(ast, StructureProxy(true, false, false, false, proxy))
    return true
end

@doc """ Post-processing function. Propagate upper part of a matrix. """
function propagate_upper_structure(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple
)
    proxy = ast.args[2]
    set_structure_proxy(ast, StructureProxy(false, true, false, false, proxy))
    return true
end

@doc """ Post-processing function. Memoize a diagonal matrix. """
function memoize_diagonal_structure(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple
)
    set_structure_proxy(ast, StructureProxy(false, false, true, false, nothing))
    return true
end

@doc """ Pre-processing function. Check if arg2 is a sparse diagonal matrix """
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
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple
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

@doc """
Pre-processing function: A function that will be called when no pattern
could handle an AST.
"""
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
    "",
    ()
)

const CS_triu_pattern = ExprPattern(
    "CS_triu_pattern",
    (:call, GlobalRef(Main, :triu), SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    propagate_upper_structure,
    "",
    "",
    ()
)

const CS_spdiagm_pattern = ExprPattern(
    "CS_spdiagm_pattern",
    (:call, GlobalRef(Main, :spdiagm), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    memoize_diagonal_structure,
    "",
    "",
    ()
)

const CS_spdiagm_times_any_pattern = ExprPattern(
    "CS_spdiagm_times_any_pattern",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    CS_spdiagm_times_any_check,
    (:NO_CHANGE, ),
    propagate_last_structure,
    "",
    "",
    ()
)

const CS_assign_pattern = ExprPattern(
    "CS_assign_pattern",
    (:(=), Any, Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    propagate_last_structure,
    "",
    "",
    ()
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
    "",                       # Useless
    ()
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


@doc """ Figure out the constant_structured property of all the matrices in the region.
"""
function constant_property_match(
    matrics     :: Dict
    region      :: LoopRegion,
    liveness    :: Liveness,
    symbol_info :: Sym2TypeMap,
    cfg         :: CFG
)
    constants   = find_constant_values(region, liveness, cfg)
    single_defs = find_single_defs(region, liveness, cfg)


    dprintln(1, 0, "\nSingle Defs:")
    dprintln(1, 1, "", single_defs)

#    call_sites  = CallSites(Set{CallSite}(), WholeFunction(), symbol_info,i
#                            Symexpr2PropertiesMap(),
#                            CS_structure_propagation_patterns,
#                            Vector{Action}(), Dict{Symexpr, Symbol}())

#    for (bb_idx, bb) in cfg.basic_blocks
#        for stmt in bb.statements
#            expr = stmt.expr
#            if typeof(expr) != Expr
#                continue
#            end
#            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
#        end
#    end
end 
