@doc """ 
Perform analyses for reordering. Write the intended transformation into actions.
"""
function reordering(
    actions     :: Vector{Action},
    func_ast    :: Expr, 
    symbol_info :: Dict{Union(Symbol,Integer), Type}, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops)
    
    dprintln(1, 0, "\nReordering phase starts...", "\nCFG:")
    dprintln(1, 1, cfg)
    
    regions = region_formation(func_ast, symbol_info, liveness, cfg, loop_info)
    for region in regions
        IAs = find_inter_dependent_arrays(region, symbol_info)
        reorder_region(region, actions, func_ast, symbol_info, liveness, cfg, loop_info)
    end

    dprintln(1, 0, "\nReordering actions to take:", actions)
    
    actions
end