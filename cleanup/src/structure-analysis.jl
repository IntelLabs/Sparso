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
Collect predefined maxtric property accoding from set_matrix_property statements
in a loop region.
"""
function find_predefined_properties(
    region      :: LoopRegion,
    cfg         :: CFG
)
    structure_proxies = Dict{Union{GenSym, Symbol}, StructureProxy}()

    for bb_idx in region.loop.members
        bb = cfg.basic_blocks[bb_idx]
        for stmt in bb.statements
            expr = stmt.expr
            if typeof(expr) != Expr || expr.head != :call
                continue
            end
            ast = expr
            m, func_name = resolve_call_names(ast.args)
            if func_name == "set_matrix_property"
                pmap = eval(ast.args[2])
                #dump(pmap)
                for (sym, p) in pmap
                    sp = StructureProxy()
                    if (p & SA_CONST_VALUED) != 0
                        sp.constant_valued = 3
                    end
                    if (p & SA_CONST_STRUCTURED) != 0
                        sp.constant_structured = 3
                    end
                    if (p & SA_SYMM_VALUED) != 0
                        sp.symmetric_valued = 3
                    end
                    if (p & SA_SYMM_STRUCTURED) != 0
                        sp.symmetric_structured = 3
                    end
                    if (p & SA_STRUCTURE_ONLY) != 0
                        # TODO: fill this property?
                    end
                    # add symbol to property map
                    structure_proxies[sym] = sp
                    dprintln(1, 1, "Predef ", sym, ": ", sp)
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
    region      :: LoopRegion,
    symbol_info :: Sym2TypeMap,
    liveness    :: Liveness,
    cfg         :: CFG
)
    dprintln(1, 0, "\n---- Matrix Structure Property Analysis -----\n")

    # symbol does not have a constant structure?
    structure_proxies = find_predefined_properties(region, cfg)

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
            p = structure_proxies[s]
            matrix_properties[s].constant_structured = (p.constant_structured>0)
            matrix_properties[s].is_symmetric = (p.symmetric_valued>0)
            matrix_properties[s].is_structure_symmetric = (p.symmetric_structured>0) 
            #matrix_properties[s].is_structure_only = structure_proxies[s].
        end
    end

    return matrix_properties
end
