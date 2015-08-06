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
       SA_REORDER_WHEN_BENEFICIAL

# Function interface
export set_options

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
include("reordering.jl")
include("context.jl")
include("show.jl")

end