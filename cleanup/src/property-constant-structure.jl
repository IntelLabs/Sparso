@doc """ Set is_constant_structured perperty for mat_property in a region """
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
        const skip_types = [GlobalRef, Int64, Float64, Bool, QuoteNode, ASCIIString]

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
                m, func_name = resolve_call_names(ast.args)
                # a quick hack for setfield! call
                if func_name == "setfield!" && ast.args[2].typ <: SparseMatrixCSC
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
        constants   = find_constant_values(region, liveness, cfg)
        single_defs = find_single_defs(region, liveness, cfg)

        if isa(region, LoopRegion)
            bb_idxs = region.loop.members
        else
            bb_idxs = keys(cfg.basic_blocks)
        end

        #sym_non_constant = gensym("SYMBOL_NON_CONSTANT")
        const sym_non_constant = Symbol(:NON_CONSTANT)

        # dependence map: k -> symbols that k depends on
        depend_map = Dict{Union{GenSym,Symbol}, Set{Union{GenSym,Symbol}}}()
        depend_map[sym_non_constant] = Set{Union{GenSym, Symbol}}()

        call_sites  = CallSites(Set{CallSite}(), region, symbol_info,
                            [],
                            Vector{Action}(), depend_map)


        # fill the dependence map by walking through all statements in the region
        for bb_idx in bb_idxs
            bb = cfg.basic_blocks[bb_idx]
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

        # property value: 
        #  -1: not constant 
        #  0: unknow 
        #  1: constant 
        #  2: external(constant)
        #  3: specified by set_matrix_property statement
        #  4: inherited from parent region (constant)
        property_map = Dict{Union{GenSym,Symbol}, Int}()
        property_map[sym_non_constant] = -1

        for k in keys(depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) && mat_property[k].constant_structured!=0 
                property_map[k] = mat_property[k].constant_structured 
            else
                property_map[k] = in(k, single_defs) ? 0 : -1
            end
        end

        for rk in keys(reverse_depend_map)
            if !haskey(property_map, rk)
                property_map[rk] = 2
            end
        end

        for c in constants
            if !haskey(property_map, c) || property_map[c] <=0
                property_map[c] = 2
            end
        end

        print_property_map(1, property_map, depend_map)

        # propagate non-constant property 
        # 
        working_set = []
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
                    if property_map[d] == 3
                        # always favor annotation
                        dprintln(1, 1, "WW annotation overwrites non_constant ", d)
                    elseif property_map[d] == 4
                        # always favor inherit result
                        dprintln(1, 1, "WW inherit overwrites non_constant ", d)
                    else
                        property_map[d] = -1
                        push!(working_set, d)
                    end
                end
            end
        end

        dprintln(1, 1, "\nAfter non-constant propagation:")
        print_property_map(1, property_map, depend_map)


        # propagate constant property 
        working_set = []
        for (k, v) in property_map
            if v > 0
                push!(working_set, k)
            end
        end

        while !isempty(working_set) 
            s = shift!(working_set)
            if !haskey(reverse_depend_map, s)
                continue
            end
            # check every symbol that depends on s
            # 
            for rd in reverse_depend_map[s]
                if property_map[rd] != 0
                    continue
                end
                constant = true
                for d in depend_map[rd]
                    if property_map[d] == 0 
                        constant = false
                    end
                end
                if constant
                    property_map[rd] = 1
                    push!(working_set, rd)
                end
            end
        end

        dprintln(1, 1, "\nAfter constant propagation:")
        print_property_map(1, property_map, depend_map)

        # copy back result
        delete!(property_map, :NON_CONSTANT)
        for (k, v) in property_map
            if any(t -> symbol_info[k]<:t, MATRIX_RELATED_TYPES)
                if !haskey(mat_property, k)
                    mat_property[k] = StructureProxy()
                end
                mat_property[k].constant_structured = v
            else
                dprintln(1, 1, "WW skip ", k, " ", symbol_info[k])
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
