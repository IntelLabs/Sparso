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

PointwiseMultiply_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseMultiply",
    (AbstractVector, AbstractVector), # The two parameters must be two vectors
    UPDATED_NONE,                     # No parameter is updated
    true,                             # The function is distributive
    ia(Set(0, 1, 2))                  # The return value (0) and the two parameters (1, 2) are inter-dependent
)

function_descriptions  = [
    PointwiseMultiply_Desc
]

println(",,,,,, here is :", function_descriptions)