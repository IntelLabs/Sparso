module SparseAccelerator

using CompilerTools
using CompilerTools.CFGs
using CompilerTools.LambdaHandling
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using CompilerTools.AstWalker
using CompilerTools.Loops

# Options user can set
export SA_ENABLE, SA_VERBOSE, SA_USE_JULIA, SA_USE_MKL, SA_USE_SPMP, 
       SA_REPLACE_CALLS, SA_REORDER, SA_CONTEXT

export @acc

# Function interface
export set_options, matrix_market_read

# Data structures user can extend
export FunctionDescription

# Explicitly specify matrix property
export  SA_CONST_VALUED, SA_CONST_STUCTURED, SA_SYMMETRIC,
        SA_STRUCTURE_SYMMETRIC, SA_STRUCTURE_ONLY

export set_matrix_property, unset_matrix_property


include("globals.jl")
include("exceptions.jl")
include("function-descriptions.jl")
include("region-formation.jl")
include("call-replacement.jl")
include("code-transformation.jl")
include("constant.jl")
include("lib-interface.jl")
include("general-reordering.jl")
include("structure-analysis.jl")
include("context.jl")
include("show.jl")

end
