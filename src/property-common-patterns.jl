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
    set_property_final_val(call_sites, ast, get_property_val(call_sites, A).final_val)
    return true
end

function prop_propagate_first_arg(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    A = ast.args[2]
    if isa(A, SymbolNode)
        A = A.name
    end
    set_property_final_val(call_sites, ast,  A)
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

    # skip if ast already has a negative property
    return vLHS.final_val != PROP_NEGATIVE_VAL
end

function prop_set_output_action(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    set_property_final_val(call_sites, ast, PROP_POSITIVE_VAL)
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

    if isempty(vLHS.vals)
        set_property_final_val(call_sites, LHS, vRHS.final_val)
    else
        assert(vLHS.final_val != PROP_NEGATIVE_VAL)

        if vRHS.final_val != vLHS.final_val
            if vRHS.final_val == PROP_NEGATIVE_VAL
                # vLHS: depend or positive, vRHS: negative
                set_property_final_val(call_sites, LHS, PROP_NEGATIVE_VAL)
            elseif vRHS.final_val == PROP_POSITIVE_VAL
                # vLHS: depend_val, vRHS: positive -> no change
                # set_property_val(call_sites, LHS, vRHS)
            else
                # vLHS: pos, vRHS: depend_val or nothing -> back to default
                set_property_final_val(call_sites, LHS, vRHS.final_val)
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
