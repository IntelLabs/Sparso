@doc """ 
Discover context information. Write the intended transformation into actions.
"""
function context_info_discovery(
    actions     :: Vector{Action},
    func_ast    :: Expr, 
    symbol_info :: Dict{Union(Symbol,Integer), Type}, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops)
    
    actions
end

