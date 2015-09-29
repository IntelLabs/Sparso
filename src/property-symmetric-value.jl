@doc """ Set is_constant_structured perperty for mat_property in a region """
type SymmetricValueProperty <: MatrixProperty 

    typealias PropertyMap Dict{Union{GenSym,Symbol, Expr}, Int}()

    @doc """ set_property_for method"""
    set_property_for    :: Function

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

                dump(k)
                if isa(k, Symbol)
                    dump(typeof(k))
                end
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
                        dump(ast.args[3])
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

    @doc """ Post-processing function. Propagate the structure of the last arg. """
    function propagate_last_structure(
        ast               :: Expr,
        call_sites        :: CallSites,
        fknob_creator     :: String,
        fknob_deletor     :: String,
        matrices_to_track :: Tuple
    )
        A = last(ast.args)
        structure = get_structure_proxy(call_sites, A) 
        if structure != nothing
            if ast.head == :(=)
                set_structure_proxy(call_sites, ast.args[1], structure_proxies[A])
            else
                set_structure_proxy(call_sites, ast, structure_proxies[A])
            end
        end
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
        structure = get_structure_proxy(call_sites, A) 
        if structure == nothing || !(structure.diagonal)
            return false
        end
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

    const CS_add_pattern = ExprPattern(
        "CS_add_pattern",
        (:(+), Any, Any),
        (:NO_SUB_PATTERNS,),
        CS_add_check,
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
    CS_symmetric_propagation_patterns = [
        CS_add_pattern,
        #CS_spdiagm_pattern,
        #CS_spdiagm_times_any_pattern,
        CS_assign_pattern,
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
            if haskey(depend_map, k)
                dprintln(1, level, v, "\t", k, "\t", set_to_str(depend_map[k]))
            else
                dprintln(1, level, v, "\t", k)
            end
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
        property_map = PropetyMap()
        property_map[sym_non_constant] = -1

        #single_defs = find_single_defs(region, liveness, cfg)

        for k in keys(depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].symmetric_valued!=0 
                    property_map[k] = mat_property[k].symmetric_valued
                end
            #else
            #    property_map[k] = in(k, single_defs) ? 0 : -1
            end
        end

        dprintln("\nBefore symmetric anaylsis:")
        print_property_map(1, property_map, depend_map)


        # cs_arg is attached to call_sites so that
        # pattern functions can pass back their results
        #  cs_arg[1]: property_map
        #  cs_arg[2]: whether a Sym/Gensym/Expr's property has been changed
        cs_arg = Tuple{PropertyMap, Bool}(property_map, false)

        call_sites  = CallSites(Set{CallSite}(), region, symbol_info,
                            CS_symmetric_propagation_patterns,
                            Vector{Action}(), cs_arg)

        converged = false
        while !converged
            converged = true
            for stmt in stmts
                expr = stmt.expr
                if typeof(expr) != Expr
                    continue
                end
                cs_arg[2] = false
                CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
                if cs_args[2] 
                    converged = false 
                end
            end
        end

        dprintln("\nAfter symmetric anaylsis:")
        print_property_map(1, property_map, depend_map)

        # copy back result
        delete!(property_map, :NON_CONSTANT)
        for (k, v) in property_map
            #if any(t -> symbol_info[k]<:t, MATRIX_RELATED_TYPES) 
                if !haskey(mat_property, k)
                    mat_property[k] = StructureProxy()
                end
                mat_property[k].constant_structured = v
            #else
            #    dprintln(1, 1, "WW skip ", k, " ", symbol_info[k])
            #end
        end
    end 

    @doc """ Constructor """
    function SymmetricValueProperty()
        instance = new()
        instance.set_property_for = set_property_for
        return instance 
    end
end
