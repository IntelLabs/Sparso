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

include("property-constant-structure.jl")
#include("property-symmetric-value.jl")

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
                dump(ast.args)
                assert(isa(ast.args[2], QuoteNode) && isa(ast.args[4], QuoteNode))
                sym = ast.args[2].value
                psym = ast.args[4].value
                if !haskey(structure_proxies, sym)
                    structure_proxies[sym] = StructureProxy()
                end
                if ast.args[3].name == :SA_LOWER_OF
                    structure_proxies[sym].lower_of = psym
                else
                    assert(ast.args[3].name == :SA_UPPER_OF)
                    structure_proxies[sym].upper_of = psym
                end
            else # set other properties
                pmap = eval(ast.args[2])
                #dump(pmap)
                assert(isa(pmap, Dict))
                for (sym, p) in pmap
                    sp = StructureProxy()
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
                        # TODO: fill this property?
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
    else
        region_name = "Func"
    end

    dprintln(1, 0, "\n---- Matrix Structure Property Analysis for " * region_name * " -----\n")

    matrix_properties = Symexpr2PropertiesMap()

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
            sp.lower_of = p.lower_of
            sp.upper_of = p.upper_of
        end
    else
        constants = find_constant_values(region, liveness, cfg)
        for c in constants
            if !haskey(structure_proxies, c)
                structure_proxies[c] = StructureProxy()
            end
            if structure_proxies[c] == 0
                structure_proxies[c].constant_valued = 1 
            end
        end
    end

    all_structure_properties = [
        ConstantStructureProperty()
    #    SymmetricValueProperty()
    ]

    for one_property in all_structure_properties
        one_property.set_property_for(structure_proxies, region, liveness, symbol_info, cfg)
    end

    # sort by keys 
    sorter = (x)->sort(collect(keys(x)))

    dprintln(1, 0, "\n" * region_name * " Matrix structures discovered:")
    dprintln(1, 1, sorter(structure_proxies))

    dprintln(1, 0, "\n" * region_name * " Constant structures discovered:")
    dprintln(1, 1, sorter(filter((k, v) -> v.constant_structured>0, structure_proxies)))

    dprintln(1, 0, "\n" * region_name * " Value symmetry discovered:")
    dprintln(1, 1, sorter(filter((k, v) -> v.symmetric_valued>0, structure_proxies)))

    dprintln(1, 0, "\n" * region_name * " Structure symmetry discovered:")
    dprintln(1, 1, sorter(filter((k, v) -> v.symmetric_structured>0, structure_proxies)))

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