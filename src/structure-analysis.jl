# Some patterns are for capturing structrures of matrices.

@doc """
"""
const SA_CONST_VALUED       = 1
const SA_CONST_STRUCTURED   = 2
const SA_SYMM_VALUED        = 4
const SA_SYMM_STRUCTURED    = 8
const SA_STRUCTURE_ONLY     = 16

@doc """ interface to explicitly specify matrix property in source code"""
function set_matrix_property(
    pmap    :: Dict{Symbol, Int}
) 
end

const SA_LOWER_OF = 1
const SA_UPPER_OF = 1

@doc """ specify that matrix A is the lower/upper part of B"""
function set_matrix_property(
    A       :: Symbol,
    part_of :: Int,
    B       :: Symbol
) 
end

@doc """
Describe part (lower, upper, diagonal) of a matrix, of which the structure matters.
"""
type StructureProxy
    constant_valued         :: Int
    constant_structured     :: Int
    symmetric_valued        :: Int
    symmetric_structured    :: Int
    lower_of                :: Any
    upper_of                :: Any
    proxy                   :: Any

    StructureProxy() = new(0, 0, 0, 0, nothing, nothing, nothing)
end

const MATRIX_RELATED_TYPES = [SparseMatrixCSC, SparseMatrix.CHOLMOD.Factor]

abstract MatrixProperty

# data types shared by analysis passes 
typealias PropertyMap Dict{Union{GenSym,Symbol}, Int}
#  
type ASTContextArgs
    changed         :: Bool
    property_map    :: Dict
    local_map       :: Dict
end

# symbol types unimportant to analysis
const skip_types = [GlobalRef, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString]


@doc """
"""
function build_depend_set_from_args(
    args    :: Array,
    call_sites :: CallSites,
    level      :: Int
)
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
                    dump(ast.args[3])
                    if !haskey(depend_sets, m)
                        depend_sets[m] = Set{Union{GenSym,Symbol}}()
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
This function is called by analysis passes (ASTWalker) 
to build an dependance map for symbols
"""
function build_dependence(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    return build_dependence(ast, call_sites, 1)
end


include("property-constant-structure.jl")
include("property-symmetric-value.jl")
include("property-symmetric-structure.jl")
include("property-lowerupper.jl")

@doc "print out the content of property proxies"
function dprint_property_proxies(
    pmap    :: Dict 
)
    dprintln(1, 1, "Sym : CV CS SV SS Lo Up")
    for (k, v) in pmap
        dprintln(1, 1, k, " : ", 
                    v.constant_valued, "  ", 
                    v.constant_structured, "  ",
                    v.symmetric_valued, "  ",
                    v.symmetric_structured, "  ",
                    v.lower_of, " ",
                    v.upper_of)
    end
end

@doc """
Collect predefined maxtric property accoding from set_matrix_property statements
in a loop region.
"""
function find_predefined_properties(
    region      :: Region,
    cfg         :: CFG
)
    structure_proxies = Dict{Union{GenSym, Symbol}, StructureProxy}()

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
                    structure_proxies[sym].lower_of = psym
                    dprintln(1, 1, "Predef: ", sym, " is lower of ", psym)
                else
                    assert(ast.args[3].name == :SA_UPPER_OF)
                    structure_proxies[sym].upper_of = psym
                    dprintln(1, 1, "Predef: ", sym, " is upper of ", psym)
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

@doc """ Find the properties of all the matrices in the region. 
A region is currently defined as a loop region. 
"""
function find_properties_of_matrices(
    region          :: Region,
    symbol_info     :: Sym2TypeMap,
    liveness        :: Liveness,
    cfg             :: CFG,
)
    if isa(region, LoopRegion)
        first_bb_idx = sort(collect(region.loop.members))[1]
        region_name = "Loop" * string(first_bb_idx)
        bb_idxs = region.loop.members
    else
        region_name = "Func"
        bb_idxs = keys(cfg.basic_blocks)
    end

    dprintln(1, 0, "\n---- Matrix Structure Property Analysis for " * region_name * " -----\n")

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
            if p.lower_of != nothing
                sp.lower_of = p.lower_of
            end
            if p.upper_of != nothing
                sp.upper_of = p.upper_of
            end
        end
    end

    constants = find_constant_values(region, liveness, cfg)
    for c in constants
        if !haskey(structure_proxies, c)
            structure_proxies[c] = StructureProxy()
        end
        if structure_proxies[c].constant_valued == 0
            structure_proxies[c].constant_valued = 1 
        end
    end

    # collect all statements in this region
    stmts = []
    for bb_idx in bb_idxs
        bb = cfg.basic_blocks[bb_idx]
        append!(stmts, bb.statements)
    end

    const sym_negative_property = Symbol(:NEGATIVE_PROPERTY)

    # dependence map: k -> symbols that k depends on
    depend_map = Dict{Union{GenSym,Symbol}, Set{Union{GenSym,Symbol}}}()
    depend_map[sym_negative_property] = Set{Union{GenSym, Symbol}}()

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

    dprintln(1, 0, "\nStructure proxies:")
    dprint_property_proxies(structure_proxies)

    all_structure_properties = [
        ConstantStructureProperty(),
        SymmetricValueProperty(),
        #SymmetricStructureProperty(),
        #LowerUpperProperty()
    ]

    dprintln(1, 0, "\nProperty anylsis passes:")
    for one_property in all_structure_properties
        one_property.set_property_for(structure_proxies, depend_map, region, liveness, symbol_info, cfg)
    end

    # sort by keys 
    sorter = (x) -> sort(collect(keys(x)), by=v->string(v))
    # filter out matrix type
    is_matrix_type = (x) -> any(t -> (haskey(symbol_info, x) && symbol_info[x]<:t), MATRIX_RELATED_TYPES)
    mfilter = (x) -> filter(is_matrix_type, x)

    dprintln(1, 0, "\n" * region_name * " Matrix structures discovered:")
    dprintln(1, 1, mfilter(sorter(structure_proxies)))

    dprintln(1, 0, "\n" * region_name * " Constant value discovered:")
    dprintln(1, 1, mfilter(sorter(filter((k, v) -> v.constant_valued>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Constant structures discovered:")
    dprintln(1, 1, mfilter(sorter(filter((k, v) -> v.constant_structured>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Value symmetry discovered:")
    dprintln(1, 1, mfilter(sorter(filter((k, v) -> v.symmetric_valued>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Structure symmetry discovered:")
    dprintln(1, 1, mfilter(sorter(filter((k, v) -> v.symmetric_structured>0, structure_proxies))))

    dprintln(1, 0, "\n" * region_name * " Upper/Lower matrix discovered:")
    for k in sorter(structure_proxies)
        v = structure_proxies[k]
        if v.lower_of != nothing
            dprintln(1, 1, k, " is lower of ", v.lower_of)
        end
        if v.upper_of != nothing
            dprintln(1, 1, k, " is upper of ", v.upper_of)
        end
    end

    # Merge with structure info
    single_defs = find_single_defs(region, liveness, cfg)
    symbols = union(single_defs, collect(keys(structure_proxies)))

    matrix_properties = Symexpr2PropertiesMap()

    for s in symbols
        # These symbols are not necessary related with matrices. But
        # some symbols may be of subtypes of Array, some other may be
        # not (like CHOLMOD.Factor, which is not a matrix, but is realted
        # with matrix). So we'd better not to filter out any symbol.
        matrix_properties[s] = MatrixProperties()
        if in(s, single_defs)
            matrix_properties[s].is_single_def = true
        end
        
        if haskey(structure_proxies, s) 
            p = structure_proxies[s]
            matrix_properties[s].constant_valued = (p.constant_valued>0)
            matrix_properties[s].constant_structured = (p.constant_structured>0)
            matrix_properties[s].is_symmetric = (p.symmetric_valued>0)
            matrix_properties[s].is_structure_symmetric = (p.symmetric_structured>0) 
            matrix_properties[s].lower_of = p.lower_of
            matrix_properties[s].upper_of = p.upper_of
        end
    end

    region.symbol_property = matrix_properties

    matrix_properties
end
