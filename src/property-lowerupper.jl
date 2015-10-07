@doc """ Set is_constant_structured perperty for mat_property in a region """
type LowerUpperProperty <: MatrixProperty 

    @doc "pass name"
    name                :: AbstractString

    @doc """ set_property_for method"""
    set_property_for    :: Function

    @doc """ Post-processing function. Memoize a diagonal matrix. """
    function memoize_diagonal_structure(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator     :: AbstractString,
        fknob_deletor     :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        push!(call_sites.extra.local_map[:SA_DIAGONAL], ast)
        return true
    end

    @doc """ Post-processing function for spmatmul_witheps. """
    function propagate_last_symbol_property(
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

    @doc """ Post-processing function for spmatmul_witheps. """
    function post_tril_triu_action(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        A = ast.args[2]
        set_property_val(call_sites, ast,  A)
        return true
    end


    @doc """ Post-processing function. Propagate the structure of the last arg. """
    function post_assign_action(
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

        if vRHS !=nothing && vLHS != vRHS
            set_property_val(call_sites, LHS, vRHS)
        end
        return true
    end

    @doc """ Post-processing function. Propagate the structure of the last arg. """
    function post_add_sub_action(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        args = ast.args[2:end]
        prop_vals = map(x -> get_property_val(call_sites, get_symexpr(x)), args)
        vR = get_property_val(call_sites, ast)
        if any(x -> x < 0, prop_vals) && vR >= 0 
            # TODO: check conflicts
            set_property_val(call_sites, ast, -1)
        elseif all(x -> x >0, prop_vals) && vR == 0 
            set_property_val(call_sites, ast, 1)
        end
        return true
    end

    @doc """ Post-processing function for * operation. """
    function post_multi_action(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        symbol_info = call_sites.symbol_info
        types = expr_skeleton(ast, symbol_info)[2:end]
        assert(length(ast.args)==length(types))

        type_map = Dict{Any, Any}()
        for (idx, arg) in enumerate(ast.args)
            type_map[arg] = types[idx]
        end

        # remove all scalars so that only matrics/vectors are left in args
        is_scalar_type = x -> (type_of_ast_node(x, symbol_info) <: Number)
        args = filter(x -> !is_scalar_type(x), ast.args[2:end])
        len = length(args)
        
        if len == 1
            v = get_property_val(call_sites, get_symexpr(args[1]))
            set_property_val(call_sites, ast, v) 
        elseif len == 2
        elseif len == 3
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
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        # TODO: put whatever you want as the last resort of the analysis here.
        return true
    end

    const CS_spdiagm_pattern = ExprPattern(
        "CS_spdiagm_pattern",
        (:call, GlobalRef(Main, :spdiagm), Vector),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        memoize_diagonal_structure,
        "",
        "",
        (),
        0,
        ()
    )

    @doc """ Pre-processing function. Check if arg2 is a sparse diagonal matrix """
    function CS_spdiagm_times_any_check(
        ast           :: Expr,
        call_sites    :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString
    )
        A = ast.args[2]
        dia = call_sites.extra.local_map[:SA_DIAGONAL] 
        #dump(A)
        #dump(dia)
        return in(A, dia)
    end

    const CS_spdiagm_times_any_pattern = ExprPattern(
        "CS_spdiagm_times_any_pattern",
        (:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        CS_spdiagm_times_any_check,
        (:NO_CHANGE, ),
        propagate_last_symbol_property,
        "",
        "",
        (),
        0,
        ()
    )

    const CS_apply_type_pattern = ExprPattern(
        "CS_apply_type_pattern",
        (:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Cdouble), GlobalRef(Main, :Cint)), Any),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        propagate_last_symbol_property,
        "",
        "",
        (),
        0,
        ()
    )


    const CS_tril_pattern = ExprPattern(
        "CS_tril_pattern",
        (:call, GlobalRef(Main, :tril), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_tril_triu_action,
        "",
        "",
        (),
        0,
        ()
    )

    const CS_triu_pattern = ExprPattern(
        "CS_triu_pattern",
        (:call, GlobalRef(Main, :triu), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_tril_triu_action,
        "",
        "",
        (),
        0,
        ()
    )


    const CS_multi_pattern = ExprPattern(
        "CS_multi_pattern",
        (:call, GlobalRef(Main, :*), Any, Any),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_multi_action,
        "",
        "",
        (),
        0,
        ()
    )

    const CS_multi3_pattern = ExprPattern(
        "CS_multi3_pattern",
        (:call, GlobalRef(Main, :*), Any, Any, Any),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_multi_action,
        "",
        "",
        (),
        0,
        ()
    )

    const CS_assign_pattern = ExprPattern(
        "CS_assign_pattern",
        (:(=), Any, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_assign_action,
        "",
        "",
        (),
        0,
        ()
    )

    const CS_assign2_pattern = ExprPattern(
        "CS_assign2_pattern",
        (:(=), SparseMatrixCSC, Any),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_assign_action,
        "",
        "",
        (),
        0,
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
        (),
        0,
        ()
    )

    @doc """
    Patterns used for discovering matrix structures. 
    """
    CS_lower_propagation_patterns = [
        CS_assign_pattern,
        CS_assign2_pattern,
        CS_tril_pattern,
        CS_spdiagm_pattern,
        CS_spdiagm_times_any_pattern,
        CS_apply_type_pattern,
        CS_last_resort_pattern
    ]

    CS_upper_propagation_patterns = [
        CS_assign_pattern,
        CS_assign2_pattern,
        CS_triu_pattern,
        CS_spdiagm_pattern,
        CS_spdiagm_times_any_pattern,
        CS_apply_type_pattern,
        CS_last_resort_pattern
    ]

    typealias PropertyKeyType   Symexpr

    const DEFAULT_PROP_VAL = nothing
    const NEG_PROP_VAL = :NEGATIVE_PROPERTY

    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        region_info         :: RegionInfo,
        mat_property        :: Dict,
    )
        property_upper_map = new_sym_property_map(Any, DEFAULT_PROP_VAL, NEG_PROP_VAL)
        property_lower_map = new_sym_property_map(Any, DEFAULT_PROP_VAL, NEG_PROP_VAL)

        for k in keys(region_info.depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].upper_of != DEFAULT_PROP_VAL
                    property_upper_map[k] = mat_property[k].upper_of
                end
                if mat_property[k].lower_of != DEFAULT_PROP_VAL
                    property_lower_map[k] = mat_property[k].lower_of
                end
            else
                property_upper_map[k] = in(k, region_info.single_defs) ? DEFAULT_PROP_VAL : NEG_PROP_VAL
            end
        end

        for rk in keys(region_info.reverse_depend_map)
            if !haskey(property_lower_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].lower_of != DEFAULT_PROP_VAL
                property_lower_map[rk] = mat_property[rk].lower_of
            end
            if !haskey(property_upper_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].upper_of != DEFAULT_PROP_VAL
                property_upper_map[rk] = mat_property[rk].upper_of
            end
        end

        dprintln(1, 0, "\nLower upper anaylsis:")

        for (pname, property_map, CS_propagation_patterns) in 
            [("lower_of", property_lower_map, CS_lower_propagation_patterns),
             ("upper_of", property_upper_map, CS_upper_propagation_patterns)]

            dprintln(1, 1, "\nBefore " * pname * " anaylsis:")
            dprint_property_map(1, property_map)

            cnt = propagate_property(property_map, region_info, CS_propagation_patterns, 
                Any, DEFAULT_PROP_VAL, NEG_PROP_VAL,
                (ctx_args) -> ctx_args.local_map[:SA_DIAGONAL] = []
            )
 
            dprintln(1, 1, "\nAfter " * pname * " anaylsis (", cnt, " iterations):")
            dprint_property_map(1, property_map)

            # copy back result
            # delete!(property_map, :NEGATIVE_PROPERTY)
            for (k, v) in property_map
                if !haskey(mat_property, k)
                    mat_property[k] = StructureProxy()
                end
                if pname == "lower_of"
                    mat_property[k].lower_of = v
                else
                    mat_property[k].upper_of = v
                end
            end
        end
    end 

    @doc """ Constructor """
    function LowerUpperProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "lower_upper_of"
        return instance 
    end
end
