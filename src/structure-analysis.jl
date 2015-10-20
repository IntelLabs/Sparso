# --- begin: set_matrix_property inferface --->

@doc """
"""
const SA_CONST_VALUED       = 1
const SA_CONST_STRUCTURED   = 2
const SA_SYMM_VALUED        = 4
const SA_SYMM_STRUCTURED    = 8
const SA_STRUCTURE_ONLY     = 16

@doc """interface to explicitly specify matrix property in source code"""
function set_matrix_property(
    pmap    :: Dict{Symbol, Int}
) 
end

const SA_LOWER_OF = 1
const SA_UPPER_OF = 2
const SA_TRANSPOSE_OF = 4

@doc """specify that matrix A is the lower/upper part or transpose of B"""
function set_matrix_property(
    A       :: Symbol,
    part_of :: Int,
    B       :: Symbol
) 
end

@doc """inspect matrix property at run time"""
function inspect_matrix_property(
    M       :: Symbol
) 
end

# --- end: set_matrix_property inferface --- 

const PROP_NEGATIVE_VAL = -1
const PROP_POSITIVE_VAL = 1
 
@doc """"""
type PropertyValue
    predefined  :: Bool
    vals        :: Set{Any}
    final_val   :: Any
 
    PropertyValue() = new(false, Set{Any}(), nothing)
end

@doc """
Describe part (lower, upper, diagonal) of a matrix, of which the structure matters.
"""
type MatrixPropertyValues
    constant_valued         :: PropertyValue
    constant_structured     :: PropertyValue
    symmetric_valued        :: PropertyValue
    symmetric_structured    :: PropertyValue
    structure_only          :: PropertyValue
    lower_of                :: PropertyValue
    upper_of                :: PropertyValue
    transpose_of            :: PropertyValue
    proxy                   :: Any

    MatrixPropertyValues() = new(PropertyValue(), PropertyValue(),
                            PropertyValue(), PropertyValue(),
                            PropertyValue(), PropertyValue(),
                            PropertyValue(), PropertyValue(),
                            nothing)
end

@doc "sorter for symbol indexed map"
sort_by_str_key = (x) -> sort(collect(keys(x)), by=v->string(v))

@doc "print out the content of property proxies"
function dprint_property_proxies(
    pmap    :: Dict 
)
    val_or_nothing = (x) -> x.final_val == nothing ? "-" : x.final_val
    dprintln(1, 1, "Sym : CV CS SV SS SO Lower Upper Trans")
    for k in sort_by_str_key(pmap)
        v = pmap[k]
        dprintln(1, 1, k, " : ", 
                    val_or_nothing(v.constant_valued), "  ", 
                    val_or_nothing(v.constant_structured), "  ",
                    val_or_nothing(v.symmetric_valued), "  ",
                    val_or_nothing(v.symmetric_structured), "  ",
                    val_or_nothing(v.structure_only), "  ",
                    val_or_nothing(v.lower_of), " ",
                    val_or_nothing(v.upper_of), " ",
                    val_or_nothing(v.transpose_of)
                    )
    end
end

const MATRIX_RELATED_TYPES = [SparseMatrixCSC, SparseMatrix.CHOLMOD.Factor, Factorization]

#
# data types shared by analysis passes 
#
abstract MatrixPropertyPass

@doc """
Allocate a new symbol to property map
"""
function new_sym_property_map()
    return Dict{Sym, PropertyValue}()
end

@doc """
Allocate a new symbol/expr to property map
"""
function new_symexpr_property_map()
    return Dict{Symexpr, PropertyValue}()
end

@doc """
Print out a property map.
"""
function dprint_property_map(
    level       :: Int,
    pmap        :: Dict,
)
    dprintln(1, level, "Sym \tPre Final\tVals")
    for k in sort_by_str_key(pmap)
        v = pmap[k]
        # skip default values
        if !v.predefined && v.final_val == nothing && isempty(v.vals)
            continue
        end
        dprintln(1, level, k, "\t", v.predefined, " ", v.final_val, "\t", v.vals)
    end
end

@doc """
Context used for analyzing one statement
"""
type StmtContextArgs
    changed             :: Set{Sym}              # symbols property has been changed?
    property_map        :: Dict{Sym, PropertyValue}      # symbol -> property
    local_map           :: Dict{Symexpr, PropertyValue}  # symbol|expr -> property
    pattern_match_filter :: Dict{Tuple{Expr, ExprPattern}, Int}
end

@doc """
Get property value for a symbol or expr
"""
function get_property_val(
    call_sites  :: CallSites,
    sym_name    :: Symexpr
)
    assert(isa(call_sites.extra, StmtContextArgs))
    if in(typeof(sym_name), [Expr])
        pmap = call_sites.extra.local_map 
    else
        pmap = call_sites.extra.property_map
        assert(haskey(pmap, sym_name))
    end

    if !haskey(pmap, sym_name)
        pmap[sym_name] = PropertyValue()
    end

    return pmap[sym_name]
end

@doc """
Add a value to property value set for a symbol or expr
"""
function set_property_final_val(
    call_sites  :: CallSites,
    sym_name    :: Symexpr,
    value       :: Any
)
    assert(isa(call_sites.extra, StmtContextArgs))
    if in(typeof(sym_name), [ Expr])
        pmap = call_sites.extra.local_map 
    else
        pmap = call_sites.extra.property_map
        assert(haskey(pmap, sym_name))

       #
        if pmap[sym_name].predefined 
            dprintln(1, 1, "WW: unable to change predefined value (", sym_name, " -> ", value,")")
            return 
        end

        if pmap[sym_name].final_val != value
            push!(call_sites.extra.changed, sym_name)
            push!(pmap[sym_name].vals, value)
            dprintln(1, 1, "set_property_val: ", sym_name, " -> ",  value)
        end
    end

    if !haskey(pmap, sym_name)
        pmap[sym_name] = PropertyValue()
    end

    pmap[sym_name].final_val = value
end

function is_property_negative(
    pv :: PropertyValue
)
    return pv.final_val == PROP_NEGATIVE_VAL || in(PROP_NEGATIVE_VAL, pv.vals)
end

# Symbol types unimportant to analysis
const skip_types = [GlobalRef, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString, Complex]

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
                assert(ast.args[2].typ == args_real_types[2])
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

    symbol_info         :: Sym2TypeMap
    liveness            :: Liveness
    cfg                 :: CFG

    constants           :: Set
    single_defs         :: Set

    pattern_match_filter :: Dict{Tuple{Expr, ExprPattern}, Int}

    function RegionInfo(
        region          :: Region, 
        symbol_info     :: Sym2TypeMap,
        liveness        :: Liveness,
        cfg             :: CFG,
    )
        info = new()
        info.region = region
        info.symbol_info = symbol_info
        info.liveness = liveness
        info.cfg = cfg

        info.pattern_match_filter = Dict{Tuple{Expr, ExprPattern}, Int}()

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

        call_sites  = CallSites(Set{CallSite}(), region, symbol_info,
                            [],
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

const prop_field_const_map = [
           (SA_CONST_VALUED, :constant_valued),
           (SA_CONST_STRUCTURED, :constant_structured),
           (SA_SYMM_VALUED, :symmetric_valued),
           (SA_SYMM_STRUCTURED, :symmetric_structured),
           (SA_STRUCTURE_ONLY, :structure_only) ]

@doc """
Collect predefined maxtric property accoding from set_matrix_property statements
in a loop region.
"""
function find_predefined_properties(
    region_info      :: RegionInfo,
)
    property_proxies = Dict{Sym, MatrixPropertyValues}()

    for stmt in region_info.stmts
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
                if !haskey(property_proxies, sym)
                    property_proxies[sym] = MatrixPropertyValues()
                end
                if ast.args[3].name == :SA_LOWER_OF
                    v = property_proxies[sym].lower_of
                    dprintln(1, 1, "Predef: ", sym, " is lower of ", psym)
                elseif ast.args[3].name == :SA_UPPER_OF
                    v = property_proxies[sym].upper_of
                    dprintln(1, 1, "Predef: ", sym, " is upper of ", psym)
                else
                    assert(ast.args[3].name == :SA_TRANSPOSE_OF)
                    v = property_proxies[sym].transpose_of
                    dprintln(1, 1, "Predef: ", sym, " is transpose of ", psym)
                end
                v.predefined = true
                v.final_val = get_symexpr(psym)
                push!(v.vals, v.final_val)
            else # set other properties
                pmap = eval(ast.args[2])
                #dump(pmap)
                assert(isa(pmap, Dict))
                for (sym, p) in pmap
                    if haskey(property_proxies, sym)
                        sp = property_proxies[sym]
                    else
                        sp = MatrixPropertyValues()
                    end

                    for (prop_const_val, prop_field) in prop_field_const_map 
                        if (p & prop_const_val) != 0
                            dprintln(1, 1, "Predef: ", sym, ".", prop_field)
                            sp.(prop_field).predefined = true
                            push!(sp.(prop_field).vals, PROP_POSITIVE_VAL)
                            sp.(prop_field).final_val = PROP_POSITIVE_VAL
                        end
                    end

                    #add symbol to property map
                    property_proxies[sym] = sp
                end
            end
            #replace the statement with nothing.
            stmt.expr = :(nothing)
        end
    end

    property_proxies
end


@doc """
Data-flow based property propagate
"""
function propagate_property(
    property_map            :: Dict, 
    region_info             :: RegionInfo,
    propagation_patterns    :: Array,
    ctx_arg_initializer     :: Any
)
    assert(ctx_arg_initializer == nothing || isa(ctx_arg_initializer, Function))

    # ctx_args is attached to call_sites so that
    # pattern functions can pass back their results
    ctx_args = StmtContextArgs(Set{Sym}(), 
                property_map, Dict{Expr, PropertyValue}(),
                region_info.pattern_match_filter,
                )

    call_sites  = CallSites(Set{CallSite}(), region_info.region, region_info.symbol_info,
                        propagation_patterns,
                        Vector{Action}(), ctx_args)

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
            ctx_args.local_map = new_symexpr_property_map()
            if ctx_arg_initializer != nothing
                ctx_arg_initializer(ctx_args)
            end
            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
        end
        if !isempty(ctx_args.changed)
            converged = true
            for s in ctx_args.changed 
               if property_map[s].final_val != old_pmap[s].final_val 
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
    symbol_info     :: Sym2TypeMap,
    liveness        :: Liveness,
    cfg             :: CFG,
)
    region_info = RegionInfo(region, symbol_info, liveness, cfg)
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
    property_proxies = find_predefined_properties(region_info)

    prop_fields = [:constant_valued, :constant_structured, 
                   :symmetric_valued, :symmetric_structured,
                   :structure_only, 
                   :lower_of, :upper_of, :transpose_of]

    #inherite positive properties from parent
    if isa(region, LoopRegion)
        for (sym, p) in region.parent.property_proxies
            if !haskey(property_proxies, sym)
                property_proxies[sym] = MatrixPropertyValues()
            end
            sp = property_proxies[sym]

            for prop in prop_fields
                p_val = getfield(p, prop)
                sp_val = getfield(sp, prop)
                if !sp_val.predefined && !is_property_negative(p_val)
                    sp_val.predefined = p_val.predefined
                    sp_val.final_val = p_val.final_val
                    sp_val.vals = copy(p_val.vals)
                end
            end
        end
    end

    for k in union(keys(region_info.depend_map), keys(region_info.reverse_depend_map))
        if !haskey(property_proxies, k)
            property_proxies[k] = MatrixPropertyValues()
        end
    end

    # do constant value analysis here
    for c in region_info.constants
        if !haskey(property_proxies, c)
            property_proxies[c] = MatrixPropertyValues()
        end
        if property_proxies[c].constant_valued.final_val == nothing
            property_proxies[c].constant_valued.final_val = PROP_POSITIVE_VAL
        end
    end

    dprintln(1, 0, "\nInitial structure proxies:")
    dprint_property_proxies(property_proxies)

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
        pass.set_property_for(region_info, property_proxies)
    end

    dprintln(1, 0, "\nFinal structure proxies:")
    dprint_property_proxies(property_proxies)

    # filter out matrix type
    is_matrix_type = (x) -> any(t -> (haskey(symbol_info, x) && symbol_info[x]<:t), MATRIX_RELATED_TYPES)
    mfilter = (x) -> filter(is_matrix_type, x)
    has_pos_val = (v) -> v.final_val != nothing && v.final_val != PROP_NEGATIVE_VAL

    dprintln(1, 0, "\n" * region_name * " Matrix structures discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(property_proxies)))

    dprintln(1, 0, "\n" * region_name * " Constant value discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> has_pos_val(v.constant_valued), property_proxies))))

    dprintln(1, 0, "\n" * region_name * " Constant structures discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> has_pos_val(v.constant_structured), property_proxies))))

    dprintln(1, 0, "\n" * region_name * " Value symmetry discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> has_pos_val(v.symmetric_valued), property_proxies))))

    dprintln(1, 0, "\n" * region_name * " Structure symmetry discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> has_pos_val(v.symmetric_structured), property_proxies))))

    dprintln(1, 0, "\n" * region_name * " Structure only discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> has_pos_val(v.structure_only), property_proxies))))

    dprintln(1, 0, "\n" * region_name * " Upper/Lower matrix discovered:")
    for k in sort_by_str_key(property_proxies)
        v = property_proxies[k]
        if has_pos_val(v.lower_of)
            dprintln(1, 1, k, " is lower of ", v.lower_of.final_val)
        end
        if has_pos_val(v.upper_of)
            dprintln(1, 1, k, " is upper of ", v.upper_of.final_val)
        end
    end

    dprintln(1, 0, "\n" * region_name * " Transpose matrix discovered:")
    for k in sort_by_str_key(property_proxies)
        v = property_proxies[k]
        if has_pos_val(v.transpose_of)
            dprintln(1, 1, k, " is transpose of ", v.transpose_of.final_val)
        end
    end

    # Merge with structure info
    symbols = union(region_info.single_defs, collect(keys(property_proxies)))

    matrix_properties = Symexpr2PropertiesMap()

    get_pos_val_or_bool = (v) -> isa(v.final_val, Symbol) ? v.final_val : v.final_val == PROP_POSITIVE_VAL

    for s in symbols
        ## These symbols are not necessary related with matrices. But
        ## some symbols may be of subtypes of Array, some other may be
        ## not (like CHOLMOD.Factor, which is not a matrix, but is realted
        ## with matrix). So we'd better not to filter out any symbol.
        matrix_properties[s] = MatrixProperties()
        if in(s, region_info.single_defs)
            matrix_properties[s].is_single_def = true
        end
        
        if haskey(property_proxies, s) 
            p = property_proxies[s]
            matrix_properties[s].constant_valued = get_pos_val_or_bool(p.constant_valued)
            matrix_properties[s].constant_structured = get_pos_val_or_bool(p.constant_structured)
            matrix_properties[s].is_symmetric = get_pos_val_or_bool(p.symmetric_valued)
            matrix_properties[s].is_structure_symmetric = get_pos_val_or_bool(p.symmetric_structured)
            matrix_properties[s].is_structure_only = get_pos_val_or_bool(p.structure_only)
            matrix_properties[s].lower_of = p.lower_of.final_val
            matrix_properties[s].upper_of = p.upper_of.final_val
            matrix_properties[s].transpose_of = p.transpose_of.final_val
        end
    end

    region.property_proxies = property_proxies
    region.symbol_property = matrix_properties

    matrix_properties
end
