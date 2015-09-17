# Some patterns are for capturing structrures of matrices.

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
        matrix_properties[s] = MatrixProperties(false, false, false, false, false, false)
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
    matrix_properties
end
