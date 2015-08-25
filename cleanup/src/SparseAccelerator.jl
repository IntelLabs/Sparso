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
       SA_REPLACE_CALLS, SA_REORDER, SA_REORDER_WHEN_BENEFICIAL, SA_CONTEXT

export @acc

# Function interface
export set_options, matrix_market_read

# Data structures user can extend
export FunctionDescription

include("globals.jl")
include("exceptions.jl")
include("function-descriptions.jl")
include("region-formation.jl")
include("distributivity.jl")
include("interdependent-arrays.jl")
include("reorderable-arrays.jl")
include("call-replacement.jl")
include("code-transformation.jl")
include("constant.jl")
include("reordering.jl")
include("structure.jl")
include("context.jl")
include("lib-interface.jl")
include("show.jl")

end
