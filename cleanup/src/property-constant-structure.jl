@doc """ Set constant_structured perperty for matrics in a region """
type ConstantStructureProperty <: MatrixProperty 

    @doc """ set_property_for method"""
    set_property_for    :: Function

    @doc """
    Match an expression pattern and do replacement.
    """
    function build_dependence(
        ast        :: Any,
        call_sites :: CallSites,
        level      :: Int
    )
        const skip_types = [GlobalRef, Int64, Float64, QuoteNode, Bool]

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
                for arg in ast.args[2:end]
                    arg_tp = typeof(arg)
                    if arg_tp <: Expr 
                        union!(depend_sets[k], build_dependence(arg, call_sites, level+1))
                    elseif arg_tp <: Symbol || typeof(arg) <: GenSym 
                        push!(depend_sets[k], arg)
                    elseif arg_tp <: SymbolNode  
                        push!(depend_sets[k], arg.name)
                    elseif in(arg_tp, skip_types)
                        # skip GlobalRef
                    else
                        dprintln(1, 2, typeof(arg), "\n")
                        dump(arg)
                        error("Unknown type")
                    end
                end
                dprintln(1, 1, k, " : ", depend_sets[k], "\n")
            elseif ast.head == :call 
                """ TODO: check function description to 
                1: get each arg's IO property
                2: handle dependence among args
                """
                for arg in ast.args[2:end]
                    arg_tp = typeof(arg)
                    if arg_tp <: Expr 
                        union!(ret_set, build_dependence(arg, call_sites, level+1))
                    elseif arg_tp <: Symbol || arg_tp <: GenSym 
                        push!(ret_set, arg)
                    elseif arg_tp <: SymbolNode  
                        push!(ret_set, arg.name)
                    elseif in(arg_tp, skip_types)
                        # skip GlobalRef
                    else
                        dprintln(1, 2, typeof(arg), "\n")
                        dump(arg)
                        error("Unknown type")
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

    #@doc """
    #Collapse a denpendence set so that only dependences among symbols are retained in the set.
    #"""
    #function collapse_depend_sets(
    #    depend_sets :: Dict 
    #)
    #    new_sets = depend_sets
    #    working_queue = Array{DependenceSet}[]
    #    append!(working_queue, keys(new_sets))
    #    while !isempty(working_queue) 
    #        k = shift!(working_queue)
    #        if !haskey(new_sets, k)
    #            continue
    #        end
    #
    #        if !haskey(new_sets, k) || isempty(new_sets[k].depends)
    #        end
    #    end
    #    return new_sets
    #end
    
    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        matrics     :: Dict,
        region      :: LoopRegion,
        liveness    :: Liveness,
        symbol_info :: Sym2TypeMap,
        cfg         :: CFG
    )
        constants   = find_constant_values(region, liveness, cfg)
        single_defs = find_single_defs(region, liveness, cfg)

        # dependence map: k -> symbols that k depends on
        depend_map = Dict{Union{GenSym,Symbol}, Set{Union{GenSym,Symbol}}}()

        call_sites  = CallSites(Set{CallSite}(), WholeFunction(), symbol_info,
                            Symexpr2PropertiesMap(),
                            [],
                            Vector{Action}(), Dict{Symexpr, Symbol}(), depend_map)

        # fill the dependence map by walking through all statements in the region
        for (bb_idx, bb) in cfg.basic_blocks
            for stmt in bb.statements
                expr = stmt.expr
                if typeof(expr) != Expr
                    continue
                end
                CompilerTools.AstWalker.AstWalk(expr, build_dependence, call_sites)
            end
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

        # property: 0: unknow, -1: not constant, 1: constant, 2: external(constant)
        property_map = Dict{Union{GenSym,Symbol}, Int}()
        for k in keys(depend_map)
            property_map[k] = in(k, single_defs) ? 0 : -1
        end
        for rk in keys(reverse_depend_map)
            if !haskey(property_map, rk)
                property_map[rk] = 2
            end
        end

        for c in constants
            if !haskey(property_map, c)
                error("Mismatched key")
            end

            if property_map[c] == 0
                property_map[c] = 1
            elseif property_map[c] < 0
                dprintln(1, 1, "WW constant property conflicts: ", c)
            end
        end

        for (k, v) in property_map
            if haskey(depend_map, k)
                dprintln(1, 1, v, "\t", k, "\t", depend_map[k])
            end
        end

        working_set = []
        # propagate non-constant property 
        for (k, v) in property_map
            if v < 0
                push!(working_set, k)
            end
        end

        while !isempty(working_set) 
            s = shift!(working_set)
            if !haskey(reverse_depend_map, s)
                continue
            end
            for d in reverse_depend_map[s]
                if property_map[d] >= 0
                    property_map[d] = -1
                    push!(working_set, d)
                end
            end
        end

        dprintln(1, 0, "after non-constant propagation:\n")
        for (k, v) in property_map
            if haskey(depend_map, k)
                dprintln(1, 1, v, "\t", k, "\t", depend_map[k])
            end
        end


        converged = false
        while converged == false
            converged = true
            for (bb_idx, bb) in cfg.basic_blocks
                for stmt in bb.statements
                    expr = stmt.expr
                    if typeof(expr) != Expr
                        continue
                    end
                    dprintln(1, 0, "", expr)
                    CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
                end
            end
        end
    end 

    @doc """ Constructor """
    function ConstantStructureProperty()
        instance = new()
        instance.set_property_for = set_property_for
        return instance 
    end
end
