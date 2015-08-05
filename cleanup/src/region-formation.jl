
@doc """ An exit of a region. """
type RegionExit
    from_BB       :: BasicBlock
    from_stmt_idx :: Int
    to_BB         :: BasicBlock
    to_stmt_idx   :: Int
end

@doc """ 
An interval of a region. It contains all or a part of statements in a block.
"""
type RegionInterval
    bb            :: BasicBlock
    from_stmt_idx :: Int
    to_stmt_idx   :: Int
    succs         :: Set{RegionInterval}
    preds         :: Set{RegionInterval}
    
    RegionInterval(block, from, to) = 
        new(block, from, to, Set{RegionInterval}(), Set{RegionInterval}())
end

@doc """ 
A region. Theorectically, it can be arbitrary part of the user code. In this
implementation, it has a single entry block dominating all the other blocks in
the region, and it must completely contain a loop (so that our transformations 
may be profitable).

A region is composed of intervals and exits. The intervals are connected with 
each other. Each block can belong to only 1 interval.
"""
type Region
    entry_bb        :: BasicBlock  # the entry block
    intervals       :: Vector{RegionInterval}
    exits           :: Set{RegionExit}
    bb2interval     :: Dict{BasicBlock, RegionInterval} # map from bb to the interval it is in
end

# Note: to be sure of profitability, we should require each path to pass through a loop. 
function DFSGrowRegion(first_BB :: CompilerTools.LivenessAnalysis.BasicBlock, 
                       current_BB :: CompilerTools.LivenessAnalysis.BasicBlock, 
                       start_stmt_idx :: Int64, 
                       lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                       in_loop :: Dict{Int,Bool}, 
                       visited :: Dict{Int,Bool}, 
                       has_loop :: Dict{Int,Bool}, 
                       bb_interval :: Dict{Int,RegionInterval}, 
                       exits :: Set{ExitEdge}, 
                       intervals :: Vector{RegionInterval}, 
                       symbolInfo :: Dict{Union(Symbol,Integer),Any}, 
                       loop_info :: CompilerTools.Loops.DomLoops, 
                       loop_bbs)
    if (DEBUG_LVL >= 2)
        println("\n\nDFSGrowRegion from BB ", current_BB.cfgbb.label, " stmt ", start_stmt_idx, "(", 
            1 <= start_stmt_idx && start_stmt_idx <= length(current_BB.statements) ? 
            current_BB.statements[start_stmt_idx].tls.index : "", ")");
    end
    
    # Disable the requirement of having a loop.
    # TODO: enable it in future
    loop_expected_on_every_path = false
    
    assert(!visited[current_BB.cfgbb.label])
    visited[current_BB.cfgbb.label]  = true
    has_loop[current_BB.cfgbb.label] = in_loop[current_BB.cfgbb.label]
    last_stmt_idx = length(current_BB.cfgbb.statements)
    for stmt_idx = start_stmt_idx : last_stmt_idx
        expr = current_BB.cfgbb.statements[stmt_idx].expr
        if typeof(expr) != Expr
            continue
        end
        if expr.head == :return
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, nothing
            end
            exit = ExitEdge(current_BB, stmt_idx - 1, current_BB, stmt_idx)
            push!(exits, exit)
            interval = RegionInterval(current_BB, start_stmt_idx, stmt_idx - 1)
            push!(intervals, interval)
            bb_interval[current_BB.cfgbb.label] = interval
            return true, interval
        end
        if expr.head == :throw
            throw("DFSGrowRegion: throw not handled")
            return false, nothing
        end
        distributive = checkDistributivity(expr, symbolInfo, true)
        if !distributive
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, nothing
            end
            if loop_bbs != nothing
                # The region cannot cover the whole loop
                return false, nothing
            end
            exit = ExitEdge(current_BB, stmt_idx - 1, current_BB, stmt_idx)
            push!(exits, exit)
            interval = RegionInterval(current_BB, start_stmt_idx, stmt_idx - 1)
            push!(intervals, interval)
            bb_interval[current_BB.cfgbb.label] = interval
            return true, interval
        end
    end
    interval = RegionInterval(current_BB, start_stmt_idx, last_stmt_idx) 
    push!(intervals, interval)
    bb_interval[current_BB.cfgbb.label] = interval
    
    # Current BB's statements have been scanned. Now successors
    for succ_BBcfg in current_BB.cfgbb.succs
        succ_BB = lives.basic_blocks[succ_BBcfg]
        if (DEBUG_LVL >= 2)
            println("looking at succ BB ", succ_BB.cfgbb.label, " whose dominators are ", loop_info.dom_dict[succ_BB.cfgbb.label])
        end
        
        if succ_BB.cfgbb.label == -2 || # the pseudo exit of the function
            (loop_bbs != nothing && !in(succ_BB.cfgbb.label, loop_bbs)) # out of loop TODO: allow region to include bbs out of loop. Just make sure it covers the whole loop, not part of it.
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, interval
            end
            exit = ExitEdge(current_BB, last_stmt_idx + 1, succ_BB, 0)
            push!(exits, exit)
            continue
        end            
        
        if !in(first_BB.cfgbb.label, loop_info.dom_dict[succ_BB.cfgbb.label])
            # first_BB does not dominate succ BB
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, interval
            end
            if in_loop[succ_BB.cfgbb.label]
                # we are going to insert reverse reodering on a new block between current and succ BB.
                # If succ BB is in a loop, the new block is also in a loop. We cannot affort that cost 
                return false, interval
            end
            
            exit = ExitEdge(current_BB, last_stmt_idx + 1, succ_BB, 0)
            push!(exits, exit)
            continue
        end
        
        if (visited[succ_BB.cfgbb.label])
            has_loop[current_BB.cfgbb.label] = has_loop[current_BB.cfgbb.label] || has_loop[succ_BB.cfgbb.label]
            push!(interval.succs, bb_interval[succ_BB.cfgbb.label])
            push!(bb_interval[succ_BB.cfgbb.label].preds, interval)
            continue
        end
        
        success, successor_interval = DFSGrowRegion(first_BB, succ_BB, 1, lives, in_loop, visited, has_loop, bb_interval, exits, intervals, symbolInfo, loop_info, loop_bbs)
        if !success
            return false, successor_interval
        else
            push!(interval.succs, successor_interval)
            push!(successor_interval.preds, interval)
        end
    end
    return true, interval
end

function growRegion(first_BB :: CompilerTools.CFGs.BasicBlock, 
                    mmread_stmt_idx :: Int64, 
                    lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                    in_loop :: Dict{Int,Bool}, 
                    symbolInfo :: Dict{Union(Symbol,Integer),Any}, 
                    loop_info :: CompilerTools.Loops.DomLoops, 
                    loop_bbs)
    visited  = Dict{Int, Bool}() # BB has been visited?
    has_loop = Dict{Int, Bool}() # Every path starting from the BB crosses a loop?
    bb_interval = Dict{Int, RegionInterval}()
    for (j,bb) in lives.cfg.basic_blocks 
        visited[bb.label]  = false
        has_loop[bb.label] = false
    end
    
    # Imagine there are always two invisible statements in any BB, whose statement indices are 0 and last_stmt_idx+1
    # So any real statement always has a statement before and after it. That is why we can do mmread_stmt_idx + 1 here
    # Similary, we can do any stmt_indx -1 as well.
    region = Region(lives.basic_blocks[first_BB], mmread_stmt_idx, Set{ExitEdge}(), RegionInterval[])
    success, interval = DFSGrowRegion(lives.basic_blocks[first_BB], lives.basic_blocks[first_BB], mmread_stmt_idx + 1, lives, in_loop, visited, has_loop, bb_interval, region.exits, region.intervals, symbolInfo, loop_info, loop_bbs)

    if (DEBUG_LVL >= 2)
        println("DFSGrowRegion successful?: ", success)
        println("Intervals of region found:")
        for interval in region.intervals
            show_interval(interval, 1)
            println("\tPreds: ")
            for pred in interval.preds
                show_interval(pred, 3)
            end
            println("\tSuccs: ")
            for succ in interval.succs
                show_interval(succ, 3)
            end
        end
    end
    
    if success
        return region, bb_interval
    else
        return nothing, bb_interval
    end
end

function regionFormationBasedOnLoop(L, lives, loop_info, symbolInfo)
    in_loop = BBsInLoop(lives, loop_info)
    return growRegion(lives.cfg.basic_blocks[L.head], 0, lives, in_loop, symbolInfo, loop_info, L.members)
end

@doc """ 
Form regions.
"""
function region_formation(
    func_ast    :: Expr, 
    symbol_info :: Dict{Union(Symbol,Integer), Type}, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops)

end
