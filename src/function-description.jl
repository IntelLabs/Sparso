# This file contains a description for a function, including the function's
# input/output parameters, distributivity, and inter-dependent arrays. 
# This gives compiler accurate information about the function without analyzing
# it. User can easily describe his/her own functions by adding entries to the
# vector below.

# A function is assumed to have the following form
#      f(parameter1, parameter2, ...)
# In the description, its return result is represented by number 0, 
# and its parameters by number 1, 2, and so on.
# Some parameters may be output as well, which will be described by
# the "output" field.
# ASSUMPTION: all the described is "must" information. That is, the information
# is true in all cases.
# ISSUE: depending on the specific inputs, a function "may" have different
# behavior. For example, it may update some parameters, or not update them.
# Two parameters (arrays) may be inter-dependent in this case, or not in another
# case. The function may be distributive in this case, but not in another.
# How to represent and use some "may" information?

immutable FunctionDescription
    module_name     :: String # Module of the function. It is "nothing" if it is not in a Julia module, e.g. if this is a C function 
    function_name   :: String # Name of the function
    parameter_types :: Tuple  # Tuple of the function parameters' types
    output          :: Set    # The parameters updated by the function
    distributive    :: Bool   # Is this function distributive?
    IA              :: Set    # A set. Each element itself is a set of inter-dependent parameters. E.g. 
                              #Set(Any[Set(0, 1, 2), Set(4, 5)]) means that result (0) and parameter 1 and 2 are inter-dependent, and parameter 3 and 4 are inter-dependent.
end

const UPDATED_NONE = Set()
ia(args ...) = Set(Any[args ...])

const PointwiseMultiply_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseMultiply",              # SparseAccelerator.PointwiseMultiply(x::Vector, y::Vector)
    (AbstractVector, AbstractVector), # The parameters must be vectors
    UPDATED_NONE,                     # No parameter is updated
    true,                             # The function is distributive
    ia(Set(0, 1, 2))                  # The return value (0) and the two parameters (1, 2) are inter-dependent
)

const PointwiseMultiply!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseMultiply!",             # SparseAccelerator.PointwiseMultiply!(w::Vector, x::Vector, y::Vector)
    (AbstractVector, AbstractVector, AbstractVector), # The parameters must be vectors
    Set(1),                           # Parameter 1 (w) is updated
    true,                             # The function is distributive
    ia(Set(0, 1, 2, 3))               # The return value (0) and the parameters (1, 2, 3) are inter-dependent. In fact, 0 and 1 are the same array
)

const PointwiseDivide_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseDivide",                # SparseAccelerator.PointwiseDivide(x::Vector, y::Vector)
    (AbstractVector, AbstractVector), # The parameters must be vectors
    UPDATED_NONE,                     # No parameter is updated
    true,                             # The function is distributive
    ia(Set(0, 1, 2))                  # The return value (0) and the two parameters (1, 2) are inter-dependent
)

const PointwiseDivide!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseDivide!",               # SparseAccelerator.PointwiseDivide!(w::Vector, x::Vector, y::Vector)
    (AbstractVector, AbstractVector, AbstractVector), # The parameters must be vectors
    Set(1),                           # Parameter 1 (w) is updated
    true,                             # The function is distributive
    ia(Set(0, 1, 2, 3))               # The return value (0) and the parameters (1, 2, 3) are inter-dependent. In fact, 0 and 1 are the same array
)

const SpMV_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",                           # SparseAccelerator.SpMV(A::SparseMatrixCSC, x::AbstractVector)
    (SparseMatrixCSC, AbstractVector),# The parameters must be vectors
    UPDATED_NONE,                     # No parameter is updated
    true,                             # The function is distributive
    ia(Set(0, 1, 2))                  # The return value (0) and the parameters (1, 2) are inter-dependent
)

const function_descriptions  = [
    PointwiseMultiply_Desc,
    PointwiseMultiply!_Desc,
    PointwiseDivide_Desc,
    PointwiseDivide!_Desc,
    SpMV_Desc
]

function show_function_descriptions()
    println("Function descriptions: ", function_descriptions)
end

function lookup_function_description(module_name :: String, function_name :: String)
    for desc in function_descriptions
        if module_name == desc.module_name  && function_name == desc.function_name
            return desc
        end
    end
    return nothing
end