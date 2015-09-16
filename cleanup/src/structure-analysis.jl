# Some patterns are for capturing structrures of matrices.

@doc """
Describe part (lower, upper, diagonal) of a matrix, of which the structure matters.
"""
type StructureProxy
#    diagonal :: Bool
#    symmetric:: Bool # Symmetric in structure (not necessarily in value)
    constant_valued        :: Int
    constant_structured    :: Int
    symmetric_valued       :: Int
    symmetric_structured   :: Int
    lower_of  :: Any
    upper_of  :: Any
    proxy     :: Any
end


@doc """ Find the properties of all the matrices in the region. 
A region is currently defined as a loop region. 
"""
function find_properties_of_matrices(
    region      :: LoopRegion,
    liveness    :: Liveness,
    symbol_info :: Sym2TypeMap,
    cfg         :: CFG
)

    dprintln(1, 0, "\n----Matrix Property Pass-----\n")

    # symbol does not have a constant structure?
    structure_proxies = Dict{Any, StructureProxy}()

    constant_property_match(structure_proxies, region, liveness, symbol_info, cfg)
   
    dprintln(1, 0, "\nMatrix structures discovered:")
    dprintln(1, 1, "", structure_proxies)

    # These are only to cause  structure-related tests fail until structure analysis succeeds.
    # TODO: Linxiang: please copy them to appropriate places, and
    # adding printing of meaningful data, so that regression.jl can test. 
    # See the above two dprintln as an example.
    dprintln(1, 0, "\nConstant structures discovered:")
    dprintln(1, 0, "\nValue symmetry discovered:")
    dprintln(1, 0, "\nStructure symmetry discovered:")

    # Merge with structure info
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
        
        # TODO: Todd, please fill the following fields:
        # matrix_properties[s].constant_structured = ?
        # matrix_properties[s].is_symmetric = ?
        # matrix_properties[s].is_structure_symmetric = ?
        # matrix_properties[s].is_structure_only = ?
    end
    matrix_properties
end
