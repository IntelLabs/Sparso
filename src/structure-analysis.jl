# --- begin: set_matrix_property inferface --->

@doc """ interface to explicitly specify matrix property in source code"""
function set_matrix_property(
    pmap    :: Dict{Symbol, Int}
) 
end

@doc """ specify that matrix A is the lower/upper part or transpose of B"""
function set_matrix_property(
    A       :: Symbol,
    part_of :: Int,
    B       :: Symbol
) 
end

# --- end: set_matrix_property inferface --- 

@doc """
Describe part (lower, upper, diagonal) of a matrix, of which the structure matters.
"""
type StructureProxy
    constant_valued         :: Int
    constant_structured     :: Int
    symmetric_valued        :: Int
    symmetric_structured    :: Int
    structure_only          :: Int
    lower_of                :: Any
    upper_of                :: Any
    transpose_of            :: Any
    proxy                   :: Any

    StructureProxy() = new(0, 0, 0, 0, 0, nothing, nothing, nothing, nothing)
end

@doc "sorter for symbol indexed map"
sort_by_str_key = (x) -> sort(collect(keys(x)), by=v->string(v))

@doc "print out the content of property proxies"
function dprint_property_proxies(
    pmap    :: Dict 
)
    dprintln(1, 1, "Sym : CV CS SV SS SO Lower Upper Trans")
    for k in sort_by_str_key(pmap)
        v = pmap[k]
        dprintln(1, 1, k, " : ", 
                    v.constant_valued, "  ", 
                    v.constant_structured, "  ",
                    v.symmetric_valued, "  ",
                    v.symmetric_structured, "  ",
                    v.structure_only, "  ",
                    v.lower_of, " ",
                    v.upper_of, " ",
                    v.transpose_of
                    )
    end
end

const MATRIX_RELATED_TYPES = [SparseMatrixCSC, SparseMatrix.CHOLMOD.Factor, Factorization]

#
# data types shared by analysis passes 
#
abstract MatrixProperty

@doc """Map key for default property value"""
const sym_default_id = Symbol(:DEFAULT_PROPERTY)

@doc """Map key for negative property value"""
const sym_negative_id = Symbol(:NEGATIVE_PROPERTY)

@doc """
Allocate a new symbol to property map
"""
function new_sym_property_map(
    val_type        :: Type,
    default_val     :: Any,
    negative_val    :: Any 
)
    @assert(default_val == nothing || typeof(default_val) <: val_type)
    @assert(negative_val == nothing || typeof(negative_val) <: val_type)
    pmap = Dict{Sym, val_type}()
    pmap[sym_default_id] = default_val
    pmap[sym_negative_id] = negative_val
    return pmap  
end

@doc """
Allocate a new symbol/expr to property map
"""
function new_symexpr_property_map(
    val_type        :: Type,
    default_val     :: Any,
    negative_val    :: Any 
)
    @assert(typeof(default_val) <: val_type)
    @assert(typeof(negative_val) <: val_type)
    pmap = Dict{Symexpr, val_type}()
    pmap[sym_default_id] = default_val
    pmap[sym_negative_id] = negative_val
    return pmap  
end

@doc """
Print out a property map.
"""
function dprint_property_map(
    level       :: Int,
    pmap        :: Dict,
)
    for k in sort_by_str_key(pmap)
        v = pmap[k]
        dprintln(1, level, v, "\t", k)
    end
end

@doc """
Context used for analyzing one statement
"""
type StmtContextArgs{ValTp}
    changed                 :: Set{Sym}              # symbols property has been changed?
    unchangable             :: Set{Sym}              # symbols whose property cannot be changed 
    property_map            :: Dict{Sym, ValTp}      # symbol -> property
    local_map               :: Dict{Symexpr, ValTp}  # symbol|expr -> property
    live_in_before_expr     :: Set{Sym}

    # Some patterns (those contain :t1!2, etc.) may hoist some subtrees of the
    # current statement before it. That splits one statement into more than one.
    # Such patterns should be matched at the last, because otherwise, other 
    # patterns may not be able to match what they should: they cannot find the
    # subtrees to match, which are no longer in the same statement.    
    non_splittable_patterns :: Vector{ExprPattern}
    splittable_patterns     :: Vector{ExprPattern}
end

@doc """
Get property value for a symbol or expr
"""
function get_property_val(
    call_sites  :: CallSites,
    sym_name    :: Symexpr
)
    @assert(isa(call_sites.extra, StmtContextArgs))
    if in(typeof(sym_name), [Expr])
        pmap = call_sites.extra.local_map 
    else
        pmap = call_sites.extra.property_map
    end

    if haskey(pmap, sym_name)
        return pmap[sym_name]
    else
        return pmap[sym_default_id]
    end
end

@doc """
Set property value for a symbol or expr
"""
function set_property_val(
    call_sites  :: CallSites,
    sym_name    :: Symexpr,
    value       :: Any
)
    @assert(isa(call_sites.extra, StmtContextArgs))
    if in(typeof(sym_name), [ Expr])
        pmap = call_sites.extra.local_map 
        @assert(!in(sym_name, call_sites.extra.unchangable))
    else
        pmap = call_sites.extra.property_map
        # 
        if in(sym_name, call_sites.extra.unchangable)
            dprintln(1, 1, "WW: reject property change (", sym_name, " -> ", value,") because ", sym_name, " in unchangeable set")
            return 
        end
        if !haskey(pmap, sym_name) || pmap[sym_name] != value
            push!(call_sites.extra.changed, sym_name)
        end
        dprintln(1, 1, "set_property_val: ", sym_name, " -> ",  value)
    end
 
    pmap[sym_name] = value
end

function is_property_val_set(
    call_sites  :: CallSites,
    sym_name    :: Symexpr
)
    @assert(isa(call_sites.extra, StmtContextArgs))
    if in(typeof(sym_name), [Expr])
        pmap = call_sites.extra.local_map 
    else
        pmap = call_sites.extra.property_map
    end

    return haskey(pmap, sym_name)
end


# Symbol types unimportant to analysis
const skip_types = [GlobalRef, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString]

@doc """
"""
function build_depend_set_from_args(
    args    :: Array,
    call_sites :: CallSites,
    level      :: Int
)
    dep_set = Set{Sym}()
    
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
            #dump(arg)
            error("Unknown type")
        end
    end

    dep_set
end

@doc """
Build dependence map for all symbols in a region
"""
function build_dependence(
    ast        :: Any,
    call_sites :: CallSites,
    level      :: Int
)
    symbol_info = call_sites.symbol_info
    patterns    = call_sites.patterns
    depend_sets = call_sites.extra

    ret_set = Set{Sym}()

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
                #dump(k)
                error("LHS is not symbol")
            end

            if !haskey(depend_sets, k)
                depend_sets[k] = Set{Sym}() 
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
                    #dump(ast.args[3])
                    if !haskey(depend_sets, m)
                        depend_sets[m] = Set{Sym}()
                    end
                    push!(depend_sets[m], :NEGATIVE_PROPERTY)
                end
                return ret_set 
            end

            ret_set = build_depend_set_from_args(ast.args[2:end], call_sites, level)

            # check function's output set
            # an arg depends on all other args (including itself?) if it's an ouput
            func_desc = look_for_function_description(m, func_name, args_real_types[2:end])
            if func_desc != nothing
                for out in func_desc.output
                    k = ast.args[out]
                    if !haskey(depend_sets, k)
                    #    depend_sets[k] = Set{Sym}() 
                    end
                    #union!(depend_sets[k], ret_set)
                end
            end
        elseif in(ast.head, [:gotoifnot, :return])
            # skip
        else
            #dump(ast)
            error("Unhandled expr type")
        end
    end

    return ret_set
end

@doc """
This function is called by analysis passes (ASTWalker) 
to build an dependance map for symbols
"""
function build_dependence_cb(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    build_dependence(ast, call_sites, 1)
    return nothing
end

@doc """
A regions's basic information:
"""
type RegionInfo
    region              :: Region
    name                :: AbstractString
    stmts               :: Array

    depend_map          :: Dict 
    reverse_depend_map  :: Dict 

    lambda              :: LambdaInfo
    symbol_info         :: Sym2TypeMap
    liveness            :: Liveness
    cfg                 :: CFG

    constants           :: Set
    single_defs         :: Set

    function RegionInfo(
        region          :: Region,
        lambda          :: LambdaInfo,
        symbol_info     :: Sym2TypeMap,
        liveness        :: Liveness,
        cfg             :: CFG,
    )
        info = new()
        info.region = region
        info.lambda = lambda
        info.symbol_info = symbol_info
        info.liveness = liveness
        info.cfg = cfg

        if isa(region, LoopRegion)
            first_bb_idx = sort(collect(region.loop.members))[1]
            info.name = "Loop" * string(first_bb_idx)
            bb_idxs = region.loop.members
        else
            info.name = "Func"
            bb_idxs = keys(cfg.basic_blocks)
        end

        # collect all statements in this region
        info.stmts = []
        for bb_idx in bb_idxs
            bb = cfg.basic_blocks[bb_idx]
            append!(info.stmts, bb.statements)
        end

        # dependence map: k -> symbols that k depends on
        info.depend_map = Dict{Sym, Set{Sym}}()
        info.depend_map[sym_negative_id] = Set{Union{GenSym, Symbol}}()

        call_sites  = CallSites(Set{CallSite}(), region, lambda, symbol_info,
                                liveness, [],
                                Vector{Action}(), info.depend_map)

        info.constants = find_constant_values(region, liveness, cfg)
        info.single_defs = find_single_defs(region, liveness, cfg)

        # fill the dependence map by walking through all statements in the region
        for stmt in info.stmts
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end
            CompilerTools.AstWalker.AstWalk(expr, build_dependence_cb, call_sites)
        end

        # reverse dependence map: k -> symbols that depends on k 
        info.reverse_depend_map = Dict{Sym, Set{Sym}}()

        # fill reverse dependence map
        for (k, s) in info.depend_map
            for v in s
                if !haskey(info.reverse_depend_map, v)
                    info.reverse_depend_map[v] = Set{Sym}()
                end
                push!(info.reverse_depend_map[v], k)
            end
        end

        return info
    end
end

include("property-constant-structure.jl")

include("property-common-patterns.jl")
include("property-symmetric-value.jl")
include("property-symmetric-structure.jl")
include("property-structure-only.jl")
include("property-lowerupper.jl")
include("property-transpose.jl")


@doc """
Collect predefined maxtric property accoding from set_matrix_property statements
in a loop region.
"""
function find_predefined_properties(
    region      :: Region,
    cfg         :: CFG
)
    structure_proxies = Dict{Sym, StructureProxy}()

    stmts = []

    if isa(region, LoopRegion)
        for bb_idx in region.loop.members
            bb = cfg.basic_blocks[bb_idx]
            append!(stmts, bb.statements)
        end
    elseif isa(region, FunctionRegion)
        # only check the entry block
        for (bb_idx, bb) in cfg.basic_blocks
            if bb.label == -1
                stmts = bb.statements
                break
            end
        end
    else
        error("Unknown region type")
    end 

    for stmt in stmts
        expr = stmt.expr
        if typeof(expr) != Expr || expr.head != :call
            continue
        end

        ast = expr
        m, func_name = resolve_call_names(ast.args)
        if func_name == "set_matrix_property"
            if length(ast.args) == 4 # set lower_of/upper_of
                #dump(ast.args)
                assert(isa(ast.args[2], QuoteNode) && isa(ast.args[4], QuoteNode))
                sym = ast.args[2].value
                psym = ast.args[4].value
                if !haskey(structure_proxies, sym)
                    structure_proxies[sym] = StructureProxy()
                end
                if ast.args[3].name == :SA_LOWER_OF
                    structure_proxies[sym].lower_of = get_symexpr(psym)
                    dprintln(1, 1, "Predef: ", sym, " is lower of ", psym)
                elseif ast.args[3].name == :SA_UPPER_OF
                    structure_proxies[sym].upper_of = get_symexpr(psym)
                    dprintln(1, 1, "Predef: ", sym, " is upper of ", psym)
                else
                    assert(ast.args[3].name == :SA_TRANSPOSE_OF)
                    structure_proxies[sym].transpose_of = get_symexpr(psym)
                    dprintln(1, 1, "Predef: ", sym, " is transpose of ", psym)
                end
            else # set other properties
                pmap = eval(ast.args[2])
                #dump(pmap)
                assert(isa(pmap, Dict))
                for (sym, p) in pmap
                    if haskey(structure_proxies, sym)
                        sp = structure_proxies[sym]
                    else
                        sp = StructureProxy()
                    end
                    if (p & SA_CONST_VALUED) != 0
                        dprintln(1, 1, "Predef: const_valued ", sym)
                        sp.constant_valued = 3
                    end
                    if (p & SA_CONST_STRUCTURED) != 0
                        dprintln(1, 1, "Predef: const_structured ", sym)
                        sp.constant_structured = 3
                    end
                    if (p & SA_SYMM_VALUED) != 0
                        dprintln(1, 1, "Predef: symmetric_valued ", sym)
                        sp.symmetric_valued = 3
                    end
                    if (p & SA_SYMM_STRUCTURED) != 0
                        dprintln(1, 1, "Predef: symmetric_structured ", sym)
                        sp.symmetric_structured = 3
                    end
                    if (p & SA_STRUCTURE_ONLY) != 0
                        dprintln(1, 1, "Predef: structure_only ", sym)
                        sp.structure_only = 3
                    end

                    # add symbol to property map
                    structure_proxies[sym] = sp
                    #dprintln(1, 1, "Predef ", sym, ": ", sp)
                end
            end
            # replace the statement with nothing.
            stmt.expr = :(nothing)
        end
    end

    structure_proxies
end


@doc """
Data-flow based property propagate
"""
function propagate_property(
    property_map            :: Dict, 
    region_info             :: RegionInfo,
    propagation_patterns    :: Vector, #{Pattern},
    PropertyValType         :: Type,
    default_prop_value      :: Any,
    negative_prop_value     :: Any,
    ctx_arg_initializer     :: Any
)
    @assert(typeof(default_prop_value) <: PropertyValType)
    @assert(typeof(negative_prop_value) <: PropertyValType)
    @assert(ctx_arg_initializer == nothing || isa(ctx_arg_initializer, Function))

    unchangeable_set = Set{Sym}(keys(property_map))

    # ctx_args is attached to call_sites so that
    # pattern functions can pass back their results
    ctx_args = StmtContextArgs{PropertyValType}(Set{Sym}(), unchangeable_set, 
                                                property_map, Dict{Expr, PropertyValType}(), 
                                                Set{Sym}(), Vector{Pattern}(), Vector{Pattern}())

    call_sites  = CallSites(Set{CallSite}(), region_info.region, region_info.lambda, 
                            region_info.symbol_info,
                            region_info.liveness, propagation_patterns,
                            Vector{Action}(), ctx_args)

    # All patterns are non-splittable.
    call_sites.extra.non_splittable_patterns = propagation_patterns
    
    converged = false
    cnt = 0
    while !converged
        converged = true
        old_pmap = copy(property_map)
        ctx_args.changed = Set{Sym}()
        for stmt in region_info.stmts
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end
            ctx_args.local_map = new_symexpr_property_map(PropertyValType, default_prop_value, negative_prop_value)
            if ctx_arg_initializer != nothing
                ctx_arg_initializer(ctx_args)
            end
            call_sites.extra.live_in_before_expr = LivenessAnalysis.live_in(stmt, region_info.liveness)
            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
        end
        if !isempty(ctx_args.changed)
            for s in ctx_args.changed 
                old_val = haskey(old_pmap, s) ? old_pmap[s] : old_pmap[sym_default_id]
                if property_map[s] != old_val
                    converged = false 
                    break
                end
            end
        end
        cnt = cnt + 1
    end
    return cnt
end


@doc """ Find the properties of all the matrices in the region. 
A region is currently defined as a loop region. 
"""
function find_properties_of_matrices(
    region          :: Region,
    lambda          :: LambdaInfo,
    symbol_info     :: Sym2TypeMap,
    liveness        :: Liveness,
    cfg             :: CFG,
)
    region_info = RegionInfo(region, lambda, symbol_info, liveness, cfg)
    region_name = region_info.name

    dprintln(1, 0, "\n---- Matrix Structure Property Analysis for " * region_info.name * " -----\n")

    # print dependence map 
    dprintln(1, 0, "\n", region_info.name, " dependence sets:")
    for k in sort_by_str_key(region_info.depend_map)
        dprintln(1, 1, k, "\t", set_to_str(region_info.depend_map[k]))
    end

    # print reverse dependence map 
    dprintln(1, 0, "\n", region_info.name, " reverse dependence sets:")
    for k in sort_by_str_key(region_info.reverse_depend_map)
        dprintln(1, 1, k, "\t", set_to_str(region_info.reverse_depend_map[k]))
    end

    # symbol does not have a constant structure?
    structure_proxies = find_predefined_properties(region, cfg)

    # inherite properties from parent
    if isa(region, LoopRegion)
        for (sym, p) in region.parent.symbol_property
            sp = StructureProxy()
            if p.constant_valued || 
                p.constant_structured || 
                p.is_structure_symmetric ||
                p.is_symmetric || 
                p.is_structure_only || 
                p.lower_of != nothing ||
                p.upper_of != nothing
                if !haskey(structure_proxies, sym)
                    structure_proxies[sym] = StructureProxy()
                end
                sp = structure_proxies[sym]
            end
            if p.constant_valued && sp.constant_valued == 0 
                sp.constant_valued = 4 
            end
            if p.constant_structured && sp.constant_structured == 0
                sp.constant_structured = 4 
            end
            if p.is_symmetric && sp.symmetric_valued == 0
                sp.symmetric_valued = 4
            end
            if p.is_structure_symmetric && sp.symmetric_structured == 0
                sp.symmetric_structured = 4 
            end
            if p.is_structure_only && sp.structure_only == 0
                sp.structure_only = 4 
            end
            if p.lower_of != nothing
                sp.lower_of = p.lower_of
            end
            if p.upper_of != nothing
                sp.upper_of = p.upper_of
            end
            if p.transpose_of != nothing
                sp.transpose_of = p.transpose_of
            end
        end
    end

    # do constant value analysis here
    for c in region_info.constants
        if !haskey(structure_proxies, c)
            structure_proxies[c] = StructureProxy()
        end
        if structure_proxies[c].constant_valued == 0
            structure_proxies[c].constant_valued = 1 
        end
    end


    dprintln(1, 0, "\nInitial structure proxies:")
    dprint_property_proxies(structure_proxies)

    structure_property_passes = [
        ConstantStructureProperty(),
        SymmetricValueProperty(),
        SymmetricStructureProperty(),
        StructureOnlyProperty(),
        LowerUpperProperty(),
        TransposeProperty()
    ]

    for pass in structure_property_passes
        dprintln(1, 0, "\nBegin ", pass.name, " pass:")
        pass.set_property_for(region_info, structure_proxies)
    end

    # filter out matrix type
    is_matrix_type = (x) -> any(t -> (haskey(symbol_info, x) && symbol_info[x]<:t), MATRIX_RELATED_TYPES)
    mfilter = (x) -> filter(is_matrix_type, x)

    delete!(structure_proxies, sym_default_id)
    delete!(structure_proxies, sym_negative_id)
    #for (k, v) in structure_proxies
    #    dprintln(1, 1, k, " ", symbol_info[k])
    #end

    dprintln(1, 0, "\n" * region_name * " Matrix structures discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(structure_proxies)))

    dprintln(1, 0, "\n" * region_name * " Constant value discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> v.constant_valued>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Constant structures discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> v.constant_structured>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Value symmetry discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> v.symmetric_valued>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Structure symmetry discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> v.symmetric_structured>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Structure only discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> v.structure_only>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Upper/Lower matrix discovered:")
    for k in sort_by_str_key(structure_proxies)
        v = structure_proxies[k]
        if v.lower_of != nothing
            dprintln(1, 1, k, " is lower of ", v.lower_of)
        end
        if v.upper_of != nothing
            dprintln(1, 1, k, " is upper of ", v.upper_of)
        end
    end

    dprintln(1, 0, "\n" * region_name * " Transpose matrix discovered:")
    for k in sort_by_str_key(structure_proxies)
        v = structure_proxies[k]
        if v.transpose_of != nothing
            dprintln(1, 1, k, " is transpose of ", v.transpose_of)
        end
    end

    # Merge with structure info
    symbols = union(region_info.single_defs, collect(keys(structure_proxies)))

    matrix_properties = Symexpr2PropertiesMap()

    for s in symbols
        # These symbols are not necessary related with matrices. But
        # some symbols may be of subtypes of Array, some other may be
        # not (like CHOLMOD.Factor, which is not a matrix, but is realted
        # with matrix). So we'd better not to filter out any symbol.
        matrix_properties[s] = MatrixProperties()
        if in(s, region_info.single_defs)
            matrix_properties[s].is_single_def = true
        end
        
        if haskey(structure_proxies, s) 
            p = structure_proxies[s]
            matrix_properties[s].constant_valued = (p.constant_valued>0)
            matrix_properties[s].constant_structured = (p.constant_structured>0)
            matrix_properties[s].is_symmetric = (p.symmetric_valued>0)
            matrix_properties[s].is_structure_symmetric = (p.symmetric_structured>0) 
            matrix_properties[s].is_structure_only = (p.structure_only>0)
            matrix_properties[s].lower_of = p.lower_of
            matrix_properties[s].upper_of = p.upper_of
            matrix_properties[s].transpose_of = p.transpose_of
        end
    end

    region.symbol_property = matrix_properties

    matrix_properties
end
