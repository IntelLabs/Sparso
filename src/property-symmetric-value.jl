@doc """ Set is_constant_structured perperty for mat_property in a region """
type SymmetricValueProperty <: MatrixProperty 

    @doc "pass name"
    name                :: AbstractString

    @doc """ set_property_for method"""
    set_property_for    :: Function

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
        if vRHS !=0 && vLHS != vRHS
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
            set_property_val(call_sites, ast, 1)
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
                set_property_val(call_sites, ast, 1)
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

    const CS_spmatmul_witheps_pattern = ExprPattern(
        "CS_Aspmatmul_witheps_pattern",
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

    const CS_A_mul_Bc_pattern = ExprPattern(
        "CS_A_mul_Bc_pattern",
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

    const CS_add_pattern = ExprPattern(
        "CS_add_pattern",
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

    const CS_add3_pattern = ExprPattern(
        "CS_add3_pattern",
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

    const CS_sub_pattern = ExprPattern(
        "CS_sub_pattern",
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

    const CS_sub3_pattern = ExprPattern(
        "CS_sub3_pattern",
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
        (:(=), SparseMatrixCSC, SparseMatrixCSC),
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
    CS_symmetric_propagation_patterns = [
        CS_assign_pattern,
        CS_add_pattern,
        CS_add3_pattern,
        CS_sub_pattern,
        CS_sub3_pattern,
        CS_multi_pattern,
        CS_multi3_pattern,
        CS_A_mul_Bc_pattern,
        CS_spmatmul_witheps_pattern,
        CS_last_resort_pattern
    ]

    typealias PropertyValType Int
    const DEFAULT_PROP_VAL = 0
    const NEG_PROP_VAL = -1

    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        region_info         :: RegionInfo,
        mat_property        :: Dict,
    )

        # property value: 
        #  -1: not symmetric
        #  0: unknow 
        #  1: symmetric 
        #  2: external(constant)
        #  3: specified by set_matrix_property statement
        #  4: inherited from parent region
        property_map = new_sym_property_map(PropertyValType, DEFAULT_PROP_VAL, NEG_PROP_VAL)

        for k in keys(region_info.depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].symmetric_valued!=DEFAULT_PROP_VAL 
                    property_map[k] = mat_property[k].symmetric_valued
                end
            else
                property_map[k] = in(k, region_info.single_defs) ? DEFAULT_PROP_VAL : NEG_PROP_VAL 
            end
        end

        for rk in keys(region_info.reverse_depend_map)
            if !haskey(property_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].symmetric_valued>0
                property_map[rk] = mat_property[rk].symmetric_valued
            end
        end

        dprintln(1, 1, "\nBefore symmetric anaylsis:")
        dprint_property_map(1, property_map)

        cnt = propagate_property(property_map, region_info, CS_symmetric_propagation_patterns, 
                PropertyValType, DEFAULT_PROP_VAL, NEG_PROP_VAL,
                nothing
        )
    
        dprintln(1, 1, "\nAfter symmetric anaylsis (", cnt, " iterations):")
        dprint_property_map(1, property_map)

        # copy back result
        for (k, v) in property_map
            if !haskey(mat_property, k)
                mat_property[k] = StructureProxy()
            end
            mat_property[k].symmetric_valued = v
        end
    end 

    @doc """ Constructor """
    function SymmetricValueProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "symmetric_valued"
        return instance 
    end
end
