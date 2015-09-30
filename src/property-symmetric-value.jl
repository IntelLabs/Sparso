typealias PropertyMap Dict{Union{GenSym,Symbol}, Int}

type ASTContextArgs
    changed         :: Bool
    property_map    :: PropertyMap 
    local_map       :: Dict{Any, Int}
end

@doc """ Set is_constant_structured perperty for mat_property in a region """
type SymmetricValueProperty <: MatrixProperty 

    @doc """ set_property_for method"""
    set_property_for    :: Function

    const skip_types = [GlobalRef, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString, UnitRange]

    @doc """
    """
    function build_depend_set_from_args(
        args    :: Array,
        call_sites :: CallSites,
        level      :: Int
    )
        const skip_types = [GlobalRef, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString]
        dep_set = Set{Union{GenSym,Symbol}}()
        
        for arg in args
            arg_tp = typeof(arg)
            if arg_tp <: Expr 
                union!(dep_set, build_dependence(arg, call_sites, level+1))
            elseif arg_tp <: Symbol || typeof(arg) <: GenSym 
                push!(dep_set, arg)
            elseif arg_tp <: SymbolNode  
                push!(dep_set, arg.name)
            elseif in(arg_tp, skip_types)
                # skip GlobalRef
            else
                dprintln(1, 2, typeof(arg), "\n")
                dump(arg)
                error("Unknown type")
            end
        end

        dep_set
    end


    @doc """
    Build dependence map for all symbols
    """
    function build_dependence(
        ast        :: Any,
        call_sites :: CallSites,
        level      :: Int
    )
        symbol_info = call_sites.symbol_info
        patterns    = call_sites.patterns
        depend_sets = call_sites.extra

        ret_set = Set{Union{GenSym,Symbol}}()

        if typeof(ast) <: Expr
            dprintln(1, level, "-> ", ast)
            if ast.head == :(=)
                if ast.head == :(=)
                    # must be at top level?
                    if level != 1
                        error("Non-top level assignment")
                    end
                end
                k =  ast.args[1]

                if typeof(k) != Symbol && typeof(k) != GenSym
                    dprintln(1, 2, k, "\n")
                    dump(k)
                    error("LHS is not symbol")
                end

                if !haskey(depend_sets, k)
                    depend_sets[k] = Set{Union{GenSym,Symbol}}() 
                end 
                
                union!(depend_sets[k],
                    build_depend_set_from_args(ast.args[2:end], call_sites, level))

                dprintln(1, 1, k, " : ", depend_sets[k], "\n")
            elseif ast.head == :call 
                m, func_name = resolve_call_names(ast.args)

                # this is not direct type of args
                args_real_types = expr_skeleton(ast, symbol_info)[2:end]

                # a quick hack for setfield! call
                if func_name == "setfield!" && ast.args[2].typ <: SparseMatrixCSC
                    @assert(ast.args[2].typ == args_real_types[2])
                    m = ast.args[2].name
                    if ast.args[3].value != :nzval
                        if !haskey(depend_sets, m)
                            depend_sets[m] = Set{Union{GenSym,Symbol}}()
                        end
                        push!(depend_sets[m], :NON_CONSTANT)
                    end
                    return ret_set 
                end

                ret_set = build_depend_set_from_args(ast.args[2:end], call_sites, level)

                # check function's output set
                # an arg depends on all other args (including itself?) if it's an ouput
                func_desc = look_for_function_description(m, func_name, args_real_types[2:end])
                if func_desc != nothing
                    for o in func_desc.output
                        k = ast.args[o]
                        if !haskey(depend_sets, k)
                        #    depend_sets[k] = Set{Union{GenSym,Symbol}}() 
                        end
                        #union!(depend_sets[k], ret_set)
                    end
                end
            elseif in(ast.head, [:gotoifnot, :return])
                # skip
            else
                dump(ast)
                error("Unhandled expr type")
            end
        end

        return ret_set
    end

    @doc """
    Match an expression pattern and do replacement.
    """
    function build_dependence(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
        build_dependence(ast, call_sites, 1)
    end

    @doc """
    Get the structure proxy for a symbol or expr
    returns nothing if it's not found
    """
    function get_property_val(
        call_sites  :: CallSites,
        sym_name    :: Union{GenSym, Symbol, Expr, SymbolNode}
    )
        if typeof(sym_name) == SymbolNode
            sym_name = sym_name.name
        end

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
        sym_name    :: Union{GenSym, Symbol, Expr, SymbolNode},
        value       :: Int
    )
        if typeof(sym_name) == SymbolNode
            sym_name = sym_name.name
        end

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
        RHS = ast.args[2]
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
        prop_vals = map(x -> get_property_val(call_sites, x), args)
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
        assert(length(ast.args)==length(types))

        type_map = Dict{Any, Any}()
        for (idx, arg) in enumerate(ast.args)
            type_map[arg] = types[idx]
        end

        # remove all scales so that only matrics/vectors are left in args
        is_scale_type = x -> isa(x, SymbolNode)&&in(x.typ, [Int32, Int64, Float64]) || in(x, [Int32, Int64, Float64])
        args = filter(x -> !is_scale_type(x), ast.args[2:end])
        dump(args)
        len = length(args)
        
        if len == 1
            v = get_property_val(call_sites, args[1])
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
        region      :: Region,
        liveness    :: Liveness,
        symbol_info :: Sym2TypeMap,
        cfg         :: CFG
    )
        # constant structure is a subset of constant value
        for (sym, v) in mat_property
            if v.constant_valued > 0 && v.constant_structured == 0
                v.constant_structured = v.constant_valued
            end
        end

        #constants   = find_constant_values(region, liveness, cfg)

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


        const sym_non_constant = Symbol(:NON_CONSTANT)

        # dependence map: k -> symbols that k depends on
        depend_map = Dict{Union{GenSym,Symbol}, Set{Union{GenSym,Symbol}}}()
        depend_map[sym_non_constant] = Set{Union{GenSym, Symbol}}()

        call_sites  = CallSites(Set{CallSite}(), region, symbol_info,
                            [],
                            Vector{Action}(), depend_map)


        # fill the dependence map by walking through all statements in the region
        for stmt in stmts
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end
            CompilerTools.AstWalker.AstWalk(expr, build_dependence, call_sites)
        end

        # reverse dependence map: k -> symbols that depends on k 
        reverse_depend_map = Dict{Union{GenSym,Symbol}, Set{Union{GenSym,Symbol}}}()

        # fill reverse dependence map
        for (k, s) in depend_map
            for v in s
                if !haskey(reverse_depend_map, v)
                    reverse_depend_map[v] = Set{Union{GenSym,Symbol}}()
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
        property_map[sym_non_constant] = -1

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
        end

        dprintln(1, 1, "\nAfter symmetric anaylsis:")
        print_property_map(1, property_map, depend_map)

        # copy back result
        delete!(property_map, :NON_CONSTANT)
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
