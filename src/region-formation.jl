@doc """ The whole function region. Not really used so far. """
type FunctionRegion <: Region
    func_ast        :: Expr
    symbol_property :: Symexpr2PropertiesMap
    
    FunctionRegion(_func_ast) = new(_func_ast, Symexpr2PropertiesMap())
end

@doc """ 
An exit edge of a loop from a block in the loop to another block outside the loop.
"""
type LoopExit
    from_bb :: BasicBlock
    to_bb   :: BasicBlock
end

@doc """ 
A region that is composed of a loop. Speeding up loops is the focus of 
Sparse Accelerator. A loop region may have a parent region (the immediate
enclosing outer loop region, or the function region), and may have some exits,
and immediate members (the indexes of the basic blocks that belong to the
loop, but not to any inner loop of it).  
"""
type LoopRegion <: Region
    parent            :: Region
    loop              :: Loop
    exits             :: Set{LoopExit}
    symbol_property   :: Symexpr2PropertiesMap
    immediate_members :: Set{BasicBlockIndex} 
end

@doc """ 
Form a region with the loop.
"""
function loop_region_formation(
    parent            :: Region,
    L                 :: Loop,
    cfg               :: CFG,
    immediate_members :: Set{BasicBlockIndex} 
)
    region = LoopRegion(parent, L, Set{LoopExit}(), Symexpr2PropertiesMap(), immediate_members)
    blocks = cfg.basic_blocks
    for bb_index in L.members
        bb = blocks[bb_index]
        for succ in bb.succs
            if !in(succ.label, L.members)
                exit = LoopExit(bb, succ)
                push!(region.exits, exit)
            end
        end
    end
    region
end

@doc """ 
Form regions for all the outermost loops in the control flow graph.
"""
function loop_region_formation(
    parent    :: Region,
    cfg       :: CFG, 
    loop_info :: DomLoops
)
    regions = Vector{LoopRegion}()
    for L in loop_info.loops
        is_outermost      = true
        immediate_members = copy(L.members)
        for L1 in loop_info.loops
            assert(L == L1 || L.members != L1.members)
            if L != L1 && L.members < L1.members
                is_outermost = false
                break
            end
            if L != L1 && L1.members < L.members
                setdiff!(immediate_members, L1.members)
            end            
        end
        
        if is_outermost
            region = loop_region_formation(parent , L, cfg, immediate_members)
            push!(regions, region)
        end
    end
    return regions
end

@doc """ 
Form regions.
"""
function region_formation(
    parent    :: Region,
    cfg       :: CFG, 
    loop_info :: DomLoops
)
    # So far, form only loop regions.
    # TODO: (1) Connect multiple loop regions together into a bigger region 
    #       (2) Extend a region to include non-loop code
    return loop_region_formation(parent, cfg, loop_info)
end
