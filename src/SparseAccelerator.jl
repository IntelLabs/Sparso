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

module SparseAccelerator

package_path = joinpath(dirname(@__FILE__), "..", "deps")
push!(LOAD_PATH, package_path)

using CompilerTools
using CompilerTools.CFGs
using CompilerTools.LambdaHandling
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using CompilerTools.AstWalker
using CompilerTools.Loops

# Options user can set
export SA_ENABLE, SA_VERBOSE, SA_USE_JULIA, SA_USE_MKL, SA_USE_SPMP, 
       SA_REPLACE_CALLS, SA_USE_SPLITTING_PATTERNS, SA_REORDER, SA_CONTEXT,
       SA_CONTEXT_FUNC, SA_DISABLE_COLLECTIVE_STRUCTURE_PREDICTION

export @acc

# Function interface
export set_options, matrix_market_read

# Data structures user can extend
export FunctionDescription

# Explicitly specify matrix property
export  SA_CONST_VALUED, SA_CONST_SIZED, SA_MAXIMAL_STRUCTURED, 
        SA_SYMM_VALUED, SA_SYMM_STRUCTURED,
        SA_STRUCTURE_ONLY, SA_HAS_DEDICATED_MEMORY

export SA_LOWER_OF, SA_UPPER_OF, SA_TRANSPOSE_OF

export set_matrix_property

include("globals.jl")
include("exceptions.jl")
include("function-descriptions.jl")
include("region-formation.jl")
include("call-replacement.jl")
include("code-transformation.jl")
include("constant.jl")
include("lib-interface.jl")
include("general-reordering.jl")
include("context.jl")
include("show.jl")
include("structure-analysis.jl")

end
