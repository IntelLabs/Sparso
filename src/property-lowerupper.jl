@doc """ Set is_constant_structured perperty for mat_property in a region """
type LowerUpperProperty <: MatrixProperty 

    @doc """ set_property_for method"""
    set_property_for    :: Function

    @doc """
    Get the structure proxy for a symbol or expr
    returns nothing if it's not found
    """
    function get_property_val(
        call_sites  :: CallSites,
        sym_name    :: Symexpr
    )
        if in(typeof(sym_name), [Expr])
            pmap = call_sites.extra.local_map 
        else
            pmap = call_sites.extra.property_map
        end

        if !haskey(pmap, sym_name)
            pmap[sym_name] = nothing
        end
        return pmap[sym_name]
    end

    function set_property_val(
        call_sites  :: CallSites,
        sym_name    :: Symexpr,
        value       :: Any
    )
        if in(typeof(sym_name), [ Expr])
            pmap = call_sites.extra.local_map 
        else
            pmap = call_sites.extra.property_map
            if value != nothing && (!haskey(pmap, sym_name) || pmap[sym_name] != value)
                call_sites.extra.changed = true
            end
        end
        dprintln(1, 1, "update val: ", sym_name,"=", value)
        pmap[sym_name] = value
    end


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
        dprintln(1, 1,"\nasdfasdfasdf")
        dump(last(ast.args))
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


    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        mat_property:: Dict,
        depend_map  :: Dict, 
        region      :: Region,
        liveness    :: Liveness,
        symbol_info :: Sym2TypeMap,
        cfg         :: CFG
    )
        if isa(region, LoopRegion)
            bb_idxs = region.loop.members
        else
            bb_idxs = keys(cfg.basic_blocks)
        end

        # collect all statements in this region
        stmts = []
        for bb_idx in bb_idxs
            bb = cfg.basic_blocks[bb_idx]
            append!(stmts, bb.statements)
        end

        # reverse dependence map: k -> symbols that depends on k 
        reverse_depend_map = Dict{Sym, Set{Sym}}()

        # fill reverse dependence map
        for (k, s) in depend_map
            for v in s
                if !haskey(reverse_depend_map, v)
                    reverse_depend_map[v] = Set{Sym}()
                end
                push!(reverse_depend_map[v], k)
            end
        end

        single_defs = find_single_defs(region, liveness, cfg)

        property_upper_map = Dict{Sym, Any}()
        property_upper_map[:NEGATIVE_PROPERTY] = :NEGATIVE_PROPERTY

        property_lower_map = Dict{Sym, Any}()
        property_lower_map[:NEGATIVE_PROPERTY] = :NEGATIVE_PROPERTY

        for k in keys(depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].upper_of != nothing
                    property_upper_map[k] = mat_property[k].upper_of
                end
                if mat_property[k].lower_of != nothing
                    property_lower_map[k] = mat_property[k].lower_of
                end
            else
                property_upper_map[k] = in(k, single_defs) ? nothing : :NEGATIVE_PROPERTY
            end
        end

        for rk in keys(reverse_depend_map)
            if !haskey(property_lower_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].lower_of != nothing
                property_lower_map[rk] = mat_property[rk].lower_of
            end
            if !haskey(property_upper_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].upper_of != nothing
                property_upper_map[rk] = mat_property[rk].upper_of
            end
        end

        for (pname, property_map, CS_propagation_patterns) in 
            [("lower_of", property_lower_map, CS_lower_propagation_patterns),
             ("upper_of", property_upper_map, CS_upper_propagation_patterns)]

            dprintln(1, 1, "\nBefore " * pname * " anaylsis:")
            print_property_map(1, property_map, depend_map)

            # ctx_args is attached to call_sites so that
            # pattern functions can pass back their results
            ctx_args = ASTContextArgs(false, property_map, Dict{Expr, Any}())

            call_sites  = CallSites(Set{CallSite}(), region, symbol_info,
                            CS_propagation_patterns,
                            Vector{Action}(), ctx_args)

            converged = false
            cnt = 0
            while !converged
                converged = true
                for stmt in stmts
                    expr = stmt.expr
                    if typeof(expr) != Expr
                        continue
                    end
                    ctx_args.changed = false
                    ctx_args.local_map = Dict{Any, Any}()
                    ctx_args.local_map[:SA_DIAGONAL] = []
                    CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
                    if ctx_args.changed 
                        converged = false 
                    end
                end
                cnt = cnt + 1
            end

            dprintln(1, 1, "\n", cnt, " iterations.\nAfter " * pname * " anaylsis:")
            print_property_map(1, property_map, depend_map)

            # copy back result
            delete!(property_map, :NEGATIVE_PROPERTY)
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
        return instance 
    end
end
