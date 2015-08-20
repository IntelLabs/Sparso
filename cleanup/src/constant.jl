# This file contains analysis of constant (values and structures).
# So far, we assume there is no alias between different symbols.
# TODO: when alias is considered, we need add a post-processing that
# is a contamination process, to remove potential variants, due to aliases.
# But there is no need to change the basic analysis.

@doc """ Find symbols/GenSyms whose values are constant in the region. """
function find_constant_values(
    region   :: LoopRegion,
    liveness :: Liveness, 
    cfg      :: CFG
)
    constants = Set{Sym}()
    L         = region.loop
    blocks    = cfg.basic_blocks
    
    # Add all the variables that are read in the region
    for bb_idx in L.members
        bb = blocks[bb_idx]
        for stmt in bb.statements
            use = LivenessAnalysis.use(stmt, liveness)
            union!(constants, use)
        end
    end

    # Remove all the variables that are written in the region
    for bb_idx in L.members
        bb = blocks[bb_idx]
        for stmt in bb.statements
            def = LivenessAnalysis.def(stmt, liveness)
            setdiff!(constants, def)
        end
    end

    constants
end