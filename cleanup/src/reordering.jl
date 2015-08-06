@doc """ 
Find the first arrays to reorder from the function's input parameters. So far,
 the first sparse matrix paramter is regarded as the only first array to reorder.
"""
function find_first_arrays_to_reorder(
    func_ast    :: Expr, 
    symbol_info :: Symbol2TypeMap
)
    assert(func_ast.head == :lambda)
    lambda = lambdaExprToLambdaInfo(func_ast)
    FAR = Vector{Symbol}()
    for i in lambda.input_params
        if type_of_ast_node(i, symbolInfo) <: AbstractSparseMatrix
            push!(FAR, i)
            break
        end
    end
    return FAR
end

@doc """ 
Perform analyses for reordering. Write the intended transformation into actions.
"""
function reordering(
    actions     :: Vector{Action},
    func_ast    :: Expr, 
    symbol_info :: Symbol2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops)
    
    dprintln(1, 0, "\nReordering phase starts...", "\nCFG:")
    dprintln(1, 1, cfg)
    
    regions = region_formation(cfg, loop_info)
    FAR     = find_first_arrays_to_reorder(func_ast, symbol_info)
    if isempty(FAR)
        return actions
    end
    
    for region in regions
        IAs = find_inter_dependent_arrays(region, symbol_info, FAR)
        reorder_region(region, actions, func_ast, symbol_info, liveness, cfg, loop_info, FAR)
    end

    dprintln(1, 0, "\nReordering actions to take:", actions)
    
    actions
end