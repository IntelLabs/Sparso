# This file contains analysis of constant values.
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

    dprintln(1, 0, "\nConstants discovered:")
    dprintln(1, 1, "", constants)

    constants
end

@doc """
Find symbols/GenSyms that are statically defined exactly once in the region.
They are not constants, but finding them is similar to finding constants.
"""
function find_single_defs(
    region   :: LoopRegion,
    liveness :: Liveness, 
    cfg      :: CFG
)
    single_defs = Set{Sym}()
    L           = region.loop
    blocks      = cfg.basic_blocks

    for bb_idx in L.members
        bb = blocks[bb_idx]
        for stmt in bb.statements
            def            = LivenessAnalysis.def(stmt, liveness)
            total          = union(single_defs, def)
            double_defined = intersect(single_defs, def)
            single_defs    = setdiff(total, double_defined)
        end
    end

    dprintln(1, 0, "\nSingle-defs discovered:")
    dprintln(1, 1, "", single_defs)

    single_defs
end