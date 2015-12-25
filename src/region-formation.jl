@doc """ The whole function region. Not really used so far. """
type FunctionRegion <: Region
    func_ast         :: Expr
    members          :: Set{BasicBlockIndex} # All members in the function
    entry            :: Any #BasicBlock
    exit             :: Any #BasicBlock
    symbol_property  :: Symexpr2PropertiesMap
    property_proxies :: Any
    
    FunctionRegion(_func_ast) = new(_func_ast, Set{BasicBlockIndex}(), nothing,
        nothing, Symexpr2PropertiesMap(), nothing)
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
    property_proxies  :: Any
    immediate_members :: Set{BasicBlockIndex}
end

@doc """ 
Form a region with the loop.
"""
function loop_region_formation(
    parent  :: Region,
    L       :: Loop,
    cfg     :: CFG,
    members :: Set{BasicBlockIndex} 
)
    region = LoopRegion(parent, L, Set{LoopExit}(), Symexpr2PropertiesMap(), nothing, members)
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
        members = copy(L.members)
        for L1 in loop_info.loops
            assert(L == L1 || L.members != L1.members)
            if L != L1 && L.members < L1.members
                is_outermost = false
                break
            end
            if L != L1 && L1.members < L.members
                setdiff!(members, L1.members)
            end            
        end
        
        if is_outermost
            region = loop_region_formation(parent , L, cfg, members)
            push!(regions, region)
        end
    end
    return regions
end

@doc """ 
Form regions.
"""
function map_blocks_to_loop_nesting_depths(
    cfg       :: CFG,
    loop_info :: DomLoops
)
    bb2depth = Dict{BasicBlock, Int}()
    for (bb_idx, bb) in cfg.basic_blocks
        bb2depth[bb] = 0
    end
    for L in loop_info.loops
        members = copy(L.members)
        depth   = 1
        for L1 in loop_info.loops
            if L != L1 && L.members < L1.members
                depth += 1
                setdiff!(members, L1.members)
            end            
        end
        for bb_idx in members
            bb = cfg.basic_blocks[bb_idx]
            bb2depth[bb] = depth
        end
    end
    bb2depth
end

@doc """ 
Form regions.
"""
function region_formation(
    func_region :: FunctionRegion,
    cfg         :: CFG, 
    loop_info   :: DomLoops
)
    for (bb_idx, bb) in cfg.basic_blocks
        push!(func_region.members, bb_idx)
        if bb_idx == -1
            func_region.entry = bb
        end
        if bb_idx == -2
            func_region.exit = bb
        end
    end

    # So far, form only loop regions.
    # TODO: (1) Connect multiple loop regions together into a bigger region 
    #       (2) Extend a region to include non-loop code
    regions = loop_region_formation(func_region, cfg, loop_info)
    
    # Create a map from blocks to loop nesting depth. A block outside any loop has
    # a depth of 0.
    bb2depth = map_blocks_to_loop_nesting_depths(cfg, loop_info)
    
    return regions, bb2depth
end
