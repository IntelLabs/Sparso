@doc """ Set is_constant_structured perperty for mat_property in a region """
type LowerUpperProperty <: MatrixPropertyPass 

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
        push!(call_sites.extra.local_map[:SA_DIAGONAL].vals, ast)
        return true
    end

    function post_ilu_action(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString
    )
        return true
    end

    const prop_ilu_pattern = ExprPattern(
        "prop_spdiagm_pattern",
        (:call, GlobalRef(SparseAccelerator, :ilu), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_ilu_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_spdiagm_pattern = ExprPattern(
        "prop_spdiagm_pattern",
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
    function spdiagm_times_any_check(
        ast           :: Expr,
        call_sites    :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString
    )
        A = ast.args[2]
        dia = call_sites.extra.local_map[:SA_DIAGONAL].vals 
        #dump(A)
        #dump(dia)
        return in(A, dia)
    end

    const prop_spdiagm_times_any_pattern = ExprPattern(
        "prop_spdiagm_times_any_pattern",
        (:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        spdiagm_times_any_check,
        (:NO_CHANGE, ),
        prop_propagate_last_symbol_property,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_tril_pattern = ExprPattern(
        "prop_tril_pattern",
        (:call, GlobalRef(Main, :tril), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        prop_propagate_first_arg, 
        "",
        "",
        (),
        0,
        ()
    )

    const prop_triu_pattern = ExprPattern(
        "prop_triu_pattern",
        (:call, GlobalRef(Main, :triu), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        prop_propagate_first_arg, 
        "",
        "",
        (),
        0,
        ()
    )


    @doc """
    Patterns used for discovering matrix structures. 
    """
    prop_lower_propagation_patterns = [
        prop_assign_pattern,
        prop_assign2_pattern,
        prop_tril_pattern,
        prop_spdiagm_pattern,
        prop_spdiagm_times_any_pattern,
        prop_apply_type_pattern,
        prop_last_resort_pattern
    ]

    prop_upper_propagation_patterns = [
        prop_assign_pattern,
        prop_assign2_pattern,
        prop_triu_pattern,
        prop_spdiagm_pattern,
        prop_spdiagm_times_any_pattern,
        prop_apply_type_pattern,
        prop_last_resort_pattern
    ]

    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        region_info         :: RegionInfo,
        mat_property        :: Dict,
    )
        property_upper_map = new_sym_property_map()
        property_lower_map = new_sym_property_map()

        for k in keys(mat_property)
            property_upper_map[k] = mat_property[k].upper_of
            property_lower_map[k] = mat_property[k].lower_of
        end

        dprintln(1, 0, "\nLower upper analysis:")

        for (pname, property_map, prop_propagation_patterns) in 
            [("lower_of", property_lower_map, prop_lower_propagation_patterns),
             ("upper_of", property_upper_map, prop_upper_propagation_patterns)]

            dprintln(1, 1, "\nBefore " * pname * " analysis:")
            dprint_property_map(2, property_map)

            cnt = propagate_property(property_map, region_info, prop_propagation_patterns, 
                (ctx_args) -> ctx_args.local_map[:SA_DIAGONAL] = PropertyValue()
            )
 
            dprintln(1, 1, "\nAfter " * pname * " analysis (", cnt, " iterations):")
            dprint_property_map(2, property_map)
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
