# A region can be an arbitrary part of the user code.
abstract Region

@doc """ 
A region that is composed of a loop. Speeding up loops is the focus of 
Sparse Accelerator.
"""
type LoopRegion <: Region
    loop  :: Loop
end

@doc """ 
Form regions for all the outermost loops in the control flow graph.
"""
function loop_region_formation(
    cfg       :: CFG, 
    loop_info :: DomLoops
)
    regions = Vector{LoopRegion}()
    for L in loop_info.loops
        is_outermost = true
        for L1 in loop_info.loops
            assert(L == L1 || L.members != L1.members)
            if L != L && L.members < L1.members
                is_outermost = false
                break
            end
        end
        
        if is_outermost
            region = LoopRegion(L)
            push!(regions, region)
        end
    end
    return regions
end

@doc """ 
Form regions.
"""
function region_formation(
    cfg       :: CFG, 
    loop_info :: DomLoops
)
    # So far, form only loop regions.
    # TODO: (1) Connect multiple loop regions together into a bigger region 
    #       (2) Extend a region to include non-loop code
    return loop_region_formation(cfg, loop_info)
end
