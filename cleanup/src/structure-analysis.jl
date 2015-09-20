# Some patterns are for capturing structrures of matrices.

@doc """
"""
const SA_CONST_VALUED           = 1
const SA_CONST_STUCTURED        = 2
const SA_SYMMETRIC              = 4
const SA_STRUCTURE_SYMMETRIC    = 8
const SA_STRUCTURE_ONLY         = 16

@doc """ interface to explicitly specify matrix property in source code"""
function set_matrix_property(args...) 
end

function unset_matrix_property(args...) 
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

@doc """
"""
function specify_properties(
    properties  :: Dict,
    region      :: LoopRegion,
    cfg         :: CFG
)
    for bb_idx in region.loop.members
        bb = cfg.basic_blocks[bb_idx]
        for stmt in bb.statements
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            ast = expr 
            if ast.head != :call 
                continue
            end

            func = ast.args[1]
            if func.name == :(set_matrix_property) || func.name == :(unset_matrix_property)
                m = ast.args[2]
                p = properties[m]
                nv = (func.name == :(set_matrix_property))
                for arg in ast.args[3:end] 
                    if arg == SA_CONST_VALUED
                        p.constant_valued = nv
                    elseif arg == SA_CONST_STUCTURED
                        p.constant_structured = nv
                    elseif arg == SA_SYMMETRIC
                        p.is_symmetric = nv
                    elseif arg == SA_STRUCTURE_SYMMETRIC
                        p.is_structure_symmetric = nv
                    elseif arg == SA_STRUCTURE_ONLY
                        p.is_structure_only = nv
                    else
                        error("Unhanled matrix property")
                    end
                end
            end
        end
    end
end

@doc """ Find the properties of all the matrices in the region. 
A region is currently defined as a loop region. 
"""
function find_properties_of_matrices(
    region      :: LoopRegion,
    symbol_info :: Sym2TypeMap,
    liveness    :: Liveness,
    cfg         :: CFG
)

    dprintln(1, 0, "\n---- Matrix Structure Property Analysis -----\n")

    # symbol does not have a constant structure?
    structure_proxies = Dict{Union{GenSym, Symbol}, StructureProxy}()

    all_structure_properties = [
        ConstantStructureProperty()
    ]

    for one_property in all_structure_properties
        one_property.set_property_for(structure_proxies, region, liveness, symbol_info, cfg)
    end

    dprintln(1, 0, "\nMatrix structures discovered:")
    dprintln(1, 1, keys(structure_proxies))

    # These are only to cause  structure-related tests fail until structure analysis succeeds.
    dprintln(1, 0, "\nConstant structures discovered:")
    res = keys(filter((k, v) -> v.constant_structured>0, structure_proxies))
    dprintln(1, 1, res)
     
    dprintln(1, 0, "\nValue symmetry discovered:")
    dprintln(1, 0, "\nStructure symmetry discovered:")

    # Merge with structure info
    constants   = find_constant_values(region, liveness, cfg)
    single_defs = find_single_defs(region, liveness, cfg)

    symbols = union(constants, single_defs, collect(keys(structure_proxies)))
    matrix_properties = Symexpr2PropertiesMap()
    for s in symbols
        # These symbols are not necessary related with matrices. But
        # some symbols may be of subtypes of Array, some other may be
        # not (like CHOLMOD.Factor, which is not a matrix, but is realted
        # with matrix). So we'd better not to filter out any symbol.
        matrix_properties[s] = MatrixProperties()
        if in(s, constants)
            matrix_properties[s].constant_valued = true
        end
        if in(s, single_defs)
            matrix_properties[s].is_single_def = true
        end
        
        if haskey(structure_proxies, s)
            matrix_properties[s].constant_structured = (structure_proxies[s].constant_structured>0)
            #matrix_properties[s].is_symmetric = ?
            #matrix_properties[s].is_structure_symmetric = ?
            #matrix_properties[s].is_structure_only = structure_proxies[s].
        end
    end

    # explicit specification overwrite compiler analysis
    specify_properties(matrix_properties, region, cfg) 

    return matrix_properties
end
