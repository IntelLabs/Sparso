@doc """ Set is_constant_structured perperty for mat_property in a region """
type SymmetricValueProperty <: MatrixPropertyPass

    @doc "pass name"
    name                :: AbstractString

    @doc """ set_property_for method"""
    set_property_for    :: Function

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
        vR = get_property_val(call_sites, ast)
        if vR.final_val == PROP_NEGATIVE_VAL
            return true
        end

        args = ast.args[2:end]
        prop_vals = map(x -> get_property_val(call_sites, get_symexpr(x)).final_val, args)
        if any(x -> x == PROP_NEGATIVE_VAL, prop_vals)
            # TODO: check conflicts
            set_property_final_val(call_sites, ast, PROP_NEGATIVE_VAL)
        elseif all(x -> x == PROP_POSITIVE_VAL, prop_vals)
            set_property_final_val(call_sites, ast, PROP_POSITIVE_VAL)
        elseif all(x -> x != nothing, prop_vals) && vR.final_val == nothing 
            assert(length(prop_vals) == 2)
            delete!(prop_vals, PROP_POSITIVE_VAL)
            set_property_final_val(call_sites, ast, first(prop_vals))
        end
        return true
    end

    @doc """ Post-processing function. Propagate the structure of the last arg. """
    function post_A_mul_Bc_action(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        A = ast.args[2]
        B = ast.args[3]
        if isa(A, SymbolNode) && isa(B, SymbolNode) && A.name == B.name
            set_property_final_val(call_sites, ast, PROP_POSITIVE_VAL)
        end
        return true
    end

    @doc """ Post-processing function for spmatmul_witheps. """
    function post_spmatmul_witheps_action(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator :: AbstractString,
        fknob_deletor :: AbstractString,
        matrices_to_track :: Tuple,
        reordering_power  :: Int,
        reordering_FAR    :: Tuple
    )
        A = ast.args[2]
        B = ast.args[3]
        if isa(A, SymbolNode) && typeof(B) <: Expr && B.head == :call
            m, func_name = resolve_call_names(B.args)
            if func_name == "ctranspose" && B.args[2].name == A.name
                set_property_final_val(call_sites, ast, PROP_POSITIVE_VAL)
            end
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
        assert(length(ast.args) == length(types))

        type_map = Dict{Any, Any}()
        for (idx, arg) in enumerate(ast.args)
            type_map[arg] = types[idx]
        end

        # remove all scalars so that only matrics/vectors are left in args
        is_scalar_type = x -> (type_of_ast_node(x, symbol_info) <: Number)
        args = filter(x -> !is_scalar_type(x), ast.args[2:end])
        #dump(args)
        len = length(args)
        
        if len == 1
            v = get_property_val(call_sites, get_symexpr(args[1])).final_val
            set_property_final_val(call_sites, ast, v) 
        elseif len == 2
        elseif len == 3
        end

        return true
    end

    const prop_spmatmul_witheps_pattern = ExprPattern(
        "prop_Aspmatmul_witheps_pattern",
        (:call, GlobalRef(Main, :spmatmul_witheps), SparseMatrixCSC, SparseMatrixCSC, Any),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_spmatmul_witheps_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_A_mul_Bc_pattern = ExprPattern(
        "prop_A_mul_Bc_pattern",
        (:call, GlobalRef(Main, :A_mul_Bc), SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_A_mul_Bc_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_add_pattern = ExprPattern(
        "prop_add_pattern",
        (:call, GlobalRef(Main, :+), SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_add_sub_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_add3_pattern = ExprPattern(
        "prop_add3_pattern",
        (:call, GlobalRef(Main, :+), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_add_sub_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_sub_pattern = ExprPattern(
        "prop_sub_pattern",
        (:call, GlobalRef(Main, :-), SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_add_sub_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_sub3_pattern = ExprPattern(
        "prop_sub3_pattern",
        (:call, GlobalRef(Main, :-), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        post_add_sub_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_multi_pattern = ExprPattern(
        "prop_multi_pattern",
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

    const prop_multi3_pattern = ExprPattern(
        "prop_multi3_pattern",
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

    @doc """
    Patterns used for discovering matrix structures. 
    """
    prop_symmetric_propagation_patterns = [
        prop_assign_pattern,
        prop_add_pattern,
        prop_add3_pattern,
        prop_sub_pattern,
        prop_sub3_pattern,
        prop_multi_pattern,
        prop_multi3_pattern,
        prop_A_mul_Bc_pattern,
        prop_spmatmul_witheps_pattern,
        prop_last_resort_pattern
    ]

    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        region_info         :: RegionInfo,
        mat_property        :: Dict,
    )
        property_map = new_sym_property_map()

        for k in keys(mat_property)
            property_map[k] = mat_property[k].symmetric_valued
        end

        dprintln(1, 1, "\nBefore symmetric_valued analysis:")
        dprint_property_map(2, property_map)

        cnt = propagate_property(property_map, region_info, 
                    prop_symmetric_propagation_patterns, nothing)
    
        # set final values

        dprintln(1, 1, "\nAfter symmetric_valued analysis (", cnt, " iterations):")
        dprint_property_map(2, property_map)
    end 

    @doc """ Constructor """
    function SymmetricValueProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "symmetric_valued"
        return instance 
    end
end
