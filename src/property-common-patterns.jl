@doc """
Commonly used property patterns.
"""

@doc """ Post-processing function for spmatmul_witheps. """
function prop_propagate_last_symbol_property(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    A = get_symexpr(last(ast.args))
    set_property_val(call_sites, ast, get_property_val(call_sites, A))
    return true
end

@doc """ Pre-assignment check """
function prop_pre_assign_check(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString
)
    LHS = ast.args[1]
    assert(typeof(LHS) <: Sym)
    vLHS = get_property_val(call_sites, LHS)

    return vLHS != call_sites.extra.local_map[sym_negative_id]
end

function prop_set_output_action(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString
)
    set_property_val(call_sites, ast, 1)
    return true
end


@doc """ Post-processing function. Propagate the structure of the last arg. """
function prop_post_assign_action(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    LHS = ast.args[1]
    assert(typeof(LHS) <: Sym)
    RHS = get_symexpr(ast.args[2])
    vRHS = get_property_val(call_sites, RHS) 
    vLHS = get_property_val(call_sites, LHS)

    if vRHS != vLHS
        if vRHS == call_sites.extra.local_map[sym_negative_id]
            set_property_val(call_sites, LHS, vRHS)
        elseif vRHS != call_sites.extra.local_map[sym_default_id] # vRHS is positive         
            if vLHS != call_sites.extra.local_map[sym_negative_id]
                set_property_val(call_sites, LHS, vRHS)
            end
        end
    end
    return true
end


@doc """
Pre-processing function: A function that will be called when no pattern
could handle an AST.
"""
function prop_last_resort(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    # TODO: put whatever you want as the last resort of the analysis here.
    return true
end

const prop_apply_type_pattern = ExprPattern(
    "prop_apply_type_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Cdouble), GlobalRef(Main, :Cint)), Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    prop_propagate_last_symbol_property,
    "",
    "",
    (),
    0,
    ()
)

const prop_assign_pattern = ExprPattern(
    "prop_assign_pattern",
    (:(=), Any, SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    prop_pre_assign_check,
    (:NO_CHANGE, ),
    prop_post_assign_action,
    "",
    "",
    (),
    0,
    ()
)

const prop_assign2_pattern = ExprPattern(
    "prop_assign2_pattern",
    (:(=), SparseMatrixCSC, Any),
    (:NO_SUB_PATTERNS,),
    prop_pre_assign_check,
    (:NO_CHANGE, ),
    prop_post_assign_action,
    "",
    "",
    (),
    0,
    ()
)

# This is the only pattern that will always be matched, justed based on its name.
# It should always be put as the last pattern, and when it is matched, 
# the last_resort() function will be called.
const prop_last_resort_pattern = ExprPattern(
    "prop_last_resort_pattern", # Useless
    (),                       # Useless
    (:NO_SUB_PATTERNS,),      # Useless
    prop_last_resort,
    (:NO_CHANGE, ),           # Useless
    do_nothing,               # Useless
    "",                       # Useless
    "",                       # Useless
    (),
    0,
    ()
)
