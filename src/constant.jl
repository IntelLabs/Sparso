#=
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

# This file contains analysis of constant values.
# So far, we assume there is no alias between different symbols.
# TODO: when alias is considered, we need add a post-processing that
# is a contamination process, to remove potential variants, due to aliases.
# But there is no need to change the basic analysis.

@doc """ Find symbols/GenSyms whose values are constant in the region. """
function find_constant_values(
    region   :: Region,
    liveness :: Liveness, 
    cfg      :: CFG
)
    constants = Set{Sym}()
    blocks    = cfg.basic_blocks

    if isa(region, LoopRegion)
        bb_idxs = region.loop.members
    else
        bb_idxs = keys(blocks)
    end
    
    # Add all the variables that are read in the region
    for bb_idx in bb_idxs
        bb = blocks[bb_idx]
        for stmt in bb.statements
            use = LivenessAnalysis.use(stmt, liveness)
            union!(constants, use)
        end
    end

    # Remove all the variables that are written in the region
    for bb_idx in bb_idxs
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
    region   :: Region,
    liveness :: Liveness, 
    cfg      :: CFG
)
    single_defs = Set{Sym}()
    blocks      = cfg.basic_blocks

    # ISSUE: What is exactly the definition and usage of single-def?
    # Should the initial values of single_defs include the symbols live into the
    # region, even if they are not defined in the region?
    if isa(region, LoopRegion)
        bb_idxs     = region.loop.members
#        head_bb     = blocks[region.loop.head]
#        single_defs = LivenessAnalysis.live_in(head_bb, liveness)
    else
        bb_idxs     = keys(blocks)
#        entry_bb    = blocks[-1]
#        single_defs = LivenessAnalysis.live_in(entry_bb, liveness)
    end

    for bb_idx in bb_idxs
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
