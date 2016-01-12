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

@doc """inspect matrix property at run time"""
function inspect_matrix_property(
    M       :: Symbol
) 
end

# --- end: set_matrix_property inferface --- 

@doc """
A regions's basic information:
"""
type RegionInfo
    region              :: Region
    name                :: AbstractString
    stmts               :: Vector{Statement}
    bblocks             :: Dict{Int, BasicBlock}

    lambda              :: LambdaInfo
    symbol_info         :: Sym2TypeMap
    liveness            :: Liveness
    cfg                 :: CFG

    constants           :: Set
    single_defs         :: Set

    pattern_match_filter :: Dict{Tuple{Expr, ExprPattern}, Int}

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

        info.pattern_match_filter = Dict{Tuple{Expr, ExprPattern}, Int}()
        info.bblocks = Dict{Int, BasicBlock}()

        if isa(region, LoopRegion)
            first_bb_idx = sort(collect(region.loop.members))[1]
            info.name = "Loop" * string(first_bb_idx)
            bb_idxs = region.loop.members
            for i in bb_idxs
                info.bblocks[i] = cfg.basic_blocks[i]
            end
        else
            info.name = "Func"
            bb_idxs = keys(cfg.basic_blocks)
            info.bblocks = cfg.basic_blocks
        end

        # collect all statements in this region
        info.stmts = []
        for bb_idx in bb_idxs
            bb = cfg.basic_blocks[bb_idx]
            append!(info.stmts, bb.statements)
        end

        info.constants = find_constant_values(region, liveness, cfg)
        info.single_defs = find_single_defs(region, liveness, cfg)

        return info
    end
end


include("symbolic-analysis.jl")
using .SymbolicAnalysis

@doc """
Describe part (lower, upper, diagonal) of a matrix, of which the structure matters.
"""
type MatrixPropertyValues
    constant_valued         :: AbstractSymbol
#    constant_size           :: AbstractSymbol
    constant_structured     :: AbstractSymbol
    symmetric_valued        :: AbstractSymbol
    symmetric_structured    :: AbstractSymbol
    structure_only          :: AbstractSymbol
    has_dedicated_memory    :: AbstractSymbol # this array (may be a matrix or vector) has a decidated memory space, so no worry of aliases or memory allocation.
    lower_of                :: AbstractSymbol
    upper_of                :: AbstractSymbol
    transpose_of            :: AbstractSymbol
    proxy                   :: Any

    MatrixPropertyValues() = new(TOP_SYMBOL, TOP_SYMBOL,
                            TOP_SYMBOL, TOP_SYMBOL,
                            TOP_SYMBOL, TOP_SYMBOL,
                            TOP_SYMBOL, TOP_SYMBOL,
                            TOP_SYMBOL, nothing)
end

@doc "sorter for symbol indexed map"
sort_by_str_key = (x) -> sort(collect(keys(x)), by=v->string(v))

function to_string(s :: AbstractSymbol)
    if s == TOP_SYMBOL
        return "T"
    elseif s == BOTTOM_SYMBOL
        return "B"
    else
        return s.value
    end
end

@doc "print out the content of property proxies"
function dprint_property_proxies(
    pmap    :: Dict 
)
    val_or_nothing = x -> isa(x, MiddleSymbol)? x.value : "-"
    dprintln(1, 1, "Sym : CV CS SV SS SO DM Lower Upper Trans")
    for k in sort_by_str_key(pmap)
        v = pmap[k]
        dprintln(1, 1, k, " : ", 
                    to_string(v.constant_valued), "  ", 
                    to_string(v.constant_structured), "  ",
                    to_string(v.symmetric_valued), "  ",
                    to_string(v.symmetric_structured), "  ",
                    to_string(v.structure_only), "  ",
                    to_string(v.has_dedicated_memory), "  ",
                    to_string(v.lower_of), " ",
                    to_string(v.upper_of), " ",
                    to_string(v.transpose_of)
                    )
    end
end

# 
const MATRIX_RELATED_TYPES = [SparseMatrixCSC, SparseMatrix.CHOLMOD.Factor, Factorization]

# Symbol types unimportant to analysis
const SKIP_TYPES = [GlobalRef, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString, Complex]

include("structure-analysis-const-size.jl")
include("structure-analysis-transpose.jl")
include("structure-analysis-symm-value.jl")

const prop_field_const_map = [
           (SA_CONST_VALUED, :constant_valued),
           (SA_CONST_STRUCTURED, :constant_structured),
           (SA_SYMM_VALUED, :symmetric_valued),
           (SA_SYMM_STRUCTURED, :symmetric_structured),
           (SA_STRUCTURE_ONLY, :structure_only),
           (SA_HAS_DEDICATED_MEMORY, :has_dedicated_memory)]

@doc """
Collect predefined maxtric property according from set_matrix_property statements
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
                    property_proxies[sym].lower_of = MiddleSymbol(psym)
                    dprintln(1, 1, "Predef: ", sym, " is lower of ", psym)
                elseif ast.args[3].name == :SA_UPPER_OF
                    property_proxies[sym].upper_of = MiddleSymbol(psym)
                    dprintln(1, 1, "Predef: ", sym, " is upper of ", psym)
                else
                    assert(ast.args[3].name == :SA_TRANSPOSE_OF)
                    property_proxies[sym].transpose_of = MiddleSymbol(psym)
                    dprintln(1, 1, "Predef: ", sym, " is transpose of ", psym)
                end
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
                            sp.(prop_field) = MiddleSymbol(:true)
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

    # symbol does not have a constant structure?
    property_proxies = find_predefined_properties(region_info)
    predefined_symbols = keys(property_proxies)

    region.property_proxies = property_proxies

    prop_fields = [:constant_valued, :constant_structured, 
                   :symmetric_valued, :symmetric_structured,
                   :structure_only, :has_dedicated_memory,
                   :lower_of, :upper_of, :transpose_of ]

    #inherite positive properties from the parent
    if isa(region, LoopRegion)
        for (sym, p) in region.parent.property_proxies
            if !haskey(property_proxies, sym)
                property_proxies[sym] = MatrixPropertyValues()
            end
            sp = property_proxies[sym]

            for prop in prop_fields
                p_val = getfield(p, prop)
                if isa(p_val, MiddleSymbol)
                    sp.(prop) = MiddleSymbol(copy(p_val.value))
                end
            end
        end
    end

    # do constant value analysis here
    for c in region_info.constants
        if !haskey(property_proxies, c)
            property_proxies[c] = MatrixPropertyValues()
        end
        property_proxies[c].constant_valued = MiddleSymbol(:true)
    end

    dprintln(1, 0, "\nInitial structure proxies:")
    dprint_property_proxies(property_proxies)

    structure_property_passes = (
         StructureAnalysisConstSize.pass_info,
         StructureAnalysisTranspose.pass_info,
         StructureAnalysisSymmValue.pass_info,
    )

    for pass in structure_property_passes
        assert(isa(pass[1], AbstractString)) # name
        assert(isa(pass[2], Tuple)) # rules
        assert(isa(pass[3], Union{Function, Void})) # preprocess
        assert(isa(pass[4], Union{Function, Void})) # postprocess
        assert(isa(pass[5], Function)) # symbolizer

        dprintln(1, 0, "\nBegin ", pass[1], " pass:")

        rules = check_and_gen_transfer_rules(pass[2])
        analyzer = SymbolicAnalyzer(pass[1], rules, pass[5])

        # preprocess
        if pass[3] == nothing
            predefined = []
        else
            predefined = pass[3](region.property_proxies, symbol_info)
        end
        
        #
        res = run_analyzer(analyzer, region_info, predefined) 
        # creat entry for k if k is missing
        for k in keys(res)
            if !haskey(property_proxies, k)
                property_proxies[k] = MatrixPropertyValues()
            end
        end

        # postprocess
        if pass[4] != nothing
            pass[4](res, region.property_proxies, symbol_info)
        end
    end

    dprintln(1, 0, "\nFinal structure proxies:")
    dprint_property_proxies(property_proxies)


    is_matrix_type = (x) ->  any(t -> (haskey(symbol_info, x) && symbol_info[x]<:t), MATRIX_RELATED_TYPES)
    mfilter = (x) -> filter(is_matrix_type, x)
    has_pos_val = (v) -> isa(v, MiddleSymbol)

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

    dprintln(1, 0, "\n" * region_name * " Has_dedicated_memory discovered:")
    dprintln(1, 1, mfilter(sort_by_str_key(filter((k, v) -> has_pos_val(v.has_dedicated_memory), property_proxies))))

    dprintln(1, 0, "\n" * region_name * " Upper/Lower matrix discovered:")
    for k in sort_by_str_key(property_proxies)
        v = property_proxies[k]
        if has_pos_val(v.lower_of)
            dprintln(1, 1, k, " is lower of ", v.lower_of.value)
        end
        if has_pos_val(v.upper_of)
            dprintln(1, 1, k, " is upper of ", v.upper_of.value)
        end
    end

    dprintln(1, 0, "\n" * region_name * " Transpose matrix discovered:")
    for k in sort_by_str_key(property_proxies)
        v = property_proxies[k]
        if has_pos_val(v.transpose_of)
            dprintln(1, 1, k, " is transpose of ", v.transpose_of.value)
        end
    end

    # Merge with structure info
    symbols = union(region_info.single_defs, collect(keys(property_proxies)))

    matrix_properties = Symexpr2PropertiesMap()

    is_middle_symbol = (v) -> isa(v, MiddleSymbol)
    get_sym_or_nothing = (v) -> isa(v, MiddleSymbol) ? v.value : nothing

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
            matrix_properties[s].constant_valued = is_middle_symbol(p.constant_valued)
            matrix_properties[s].constant_structured = is_middle_symbol(p.constant_structured)
            matrix_properties[s].is_symmetric = is_middle_symbol(p.symmetric_valued)
            matrix_properties[s].is_structure_symmetric = is_middle_symbol(p.symmetric_structured)
            matrix_properties[s].is_structure_only = is_middle_symbol(p.structure_only)
            matrix_properties[s].has_dedicated_memory = is_middle_symbol(p.has_dedicated_memory)
            matrix_properties[s].lower_of = get_sym_or_nothing(p.lower_of)
            matrix_properties[s].upper_of = get_sym_or_nothing(p.upper_of)
            matrix_properties[s].transpose_of = get_sym_or_nothing(p.transpose_of)
        end
    end

    region.symbol_property = matrix_properties

    matrix_properties
end
