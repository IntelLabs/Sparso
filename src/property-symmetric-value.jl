@doc """ Set is_constant_structured perperty for mat_property in a region """
type SymmetricValueProperty <: MatrixProperty 

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
            pmap[sym_name] = 0
        end
        return pmap[sym_name]
    end

    function set_property_val(
        call_sites  :: CallSites,
        sym_name    :: Symexpr,
        value       :: Int
    )
        if in(typeof(sym_name), [ Expr])
            pmap = call_sites.extra.local_map 
        else
            pmap = call_sites.extra.property_map
            if value != 0 && (!haskey(pmap, sym_name) || pmap[sym_name] != value)
                call_sites.extra.changed = true
                dprintln(1, 1, "val changed: ", sym_name, value)
            end
        end
        
        pmap[sym_name] = value
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


    @doc """
    Print out a property map.
    """
    function print_property_map(
        level       :: Int,
        pm          :: Dict,
        depend_map  :: Dict
    )
        for (k, v) in pm
            #if haskey(depend_map, k)
            #    dprintln(1, level, v, "\t", k, "\t", set_to_str(depend_map[k]))
            #else
            #if isa(k, Expr)
            #else
                dprintln(1, level, v, "\t", k)
            #end
        end
    end

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

        # property value: 
        #  -1: not symmetric
        #  0: unknow 
        #  1: symmetric 
        #  2: external(constant)
        #  3: specified by set_matrix_property statement
        #  4: inherited from parent region
        property_map = PropertyMap()
        property_map[:NEGATIVE_PROPERTY] = -1

        single_defs = find_single_defs(region, liveness, cfg)

        for k in keys(depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].symmetric_valued!=0 
                    property_map[k] = mat_property[k].symmetric_valued
                end
            else
                property_map[k] = in(k, single_defs) ? 0 : -1
            end
        end

        for rk in keys(reverse_depend_map)
            if !haskey(property_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].symmetric_valued>0
                property_map[rk] = mat_property[rk].symmetric_valued
            end
        end


        dprintln(1, 1, "\nBefore symmetric anaylsis:")
        print_property_map(1, property_map, depend_map)


        # ctx_args is attached to call_sites so that
        # pattern functions can pass back their results
        ctx_args = ASTContextArgs(false, property_map, Dict{Expr, Int}())

        call_sites  = CallSites(Set{CallSite}(), region, symbol_info,
                            CS_symmetric_propagation_patterns,
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
                ctx_args.local_map = Dict{Any, Int}()
                CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
                if ctx_args.changed 
                    converged = false 
                end
            end
            cnt = cnt + 1
        end

        dprintln(1, 1, "\n", cnt, " iterations.\nAfter symmetric anaylsis:")
        print_property_map(1, property_map, depend_map)

        # copy back result
        delete!(property_map, :NEGATIVE_PROPERTY)
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
        return instance 
    end
end
