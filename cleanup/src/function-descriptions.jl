# This file contains a description for a function, including the function's
# input/output arguments, distributivity, and inter-dependent arrays. 
# This gives compiler accurate information about the function without analyzing
# it. User can easily describe his/her own functions by adding entries to the
# vector below.

# A function is assumed to have the following form
#      f(argument1, argument2, ...)
# In the description, its return result is represented by number 0, 
# and its arguments by number 1, 2, and so on.
# Some arguments may be output as well, which will be described by
# the "output" field.
# ASSUMPTION: all the described is "must" information. That is, the information
# is true in all cases.
# ISSUE: depending on the specific inputs, a function "may" have different
# behavior. For example, it may update some arguments, or not update them.
# Two arguments (arrays) may be inter-dependent in this case, or not in another
# case. The function may be distributive in this case, but not in another.
# How to represent and use some "may" information?

immutable FunctionDescription
    module_name     :: String # Module of the function. It is "nothing" if it is not in a Julia module, e.g. if this is a C function 
    function_name   :: String # Name of the function
    argument_types  :: Tuple  # Tuple of the function arguments' types
    output          :: Set    # The arguments updated by the function
    distributive    :: Bool   # Is this function distributive?
    IA              :: Set    # A set. Each element itself is a set of inter-dependent arguments. E.g. 
                              #Set(Any[Set(0, 1, 2), Set(4, 5)]) means that result (0) and argument 1 and 2 are inter-dependent, and argument 3 and 4 are inter-dependent.
end

const UPDATED_NONE  = Set()
const IA_NONE       = Set()

const VECTOR_OR_NUM = Union(Vector, Number)
const MATRIX_OR_NUM = Union(AbstractMatrix, Number)

ia(args ...) = Set(Any[args ...])

const PointwiseMultiply_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseMultiply",              # SparseAccelerator.PointwiseMultiply(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the two arguments (1, 2) are inter-dependent
)

const PointwiseMultiply!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseMultiply!",             # SparseAccelerator.PointwiseMultiply!(w::Vector, x::Vector, y::Vector)
    (Vector, Vector, Vector),         # The arguments must be vectors
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2; 3]))             # The return value (0) and the arguments (1, 2, 3) are inter-dependent. In fact, 0 and 1 are the same array
)

const PointwiseDivide_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseDivide",                # SparseAccelerator.PointwiseDivide(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the two arguments (1, 2) are inter-dependent
)

const PointwiseDivide1_Desc = FunctionDescription(
    "Main", 
    "./",           
    (VECTOR_OR_NUM, Vector),          # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the two arguments (1, 2) are inter-dependent
)

const PointwiseDivide!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "PointwiseDivide!",               # SparseAccelerator.PointwiseDivide!(w::Vector, x::Vector, y::Vector)
    (Vector, Vector, Vector),         # The arguments must be vectors
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2; 3]))             # The return value (0) and the arguments (1, 2, 3) are inter-dependent. In fact, 0 and 1 are the same array
)

const SpMV_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",                           # SparseAccelerator.SpMV(A::SparseMatrixCSC, x::Vector)
    (SparseMatrixCSC, Vector),        # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const star_Desc = FunctionDescription(
    "Main", 
    "*",                              # *(A::SparseMatrixCSC, x::Vector)
    (SparseMatrixCSC, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const star1_Desc = FunctionDescription(
    "Main", 
    "*",                              
    (Number, SparseMatrixCSC, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 2; 3]))                # The return value (0) and the arguments (2, 3) are inter-dependent
)

const star2_Desc = FunctionDescription(
    "Main", 
    "*",           
    (Number, Vector),                 
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 2]))                   # The return value (0) and the argument 2 are inter-dependent
)

const Dot_Desc = FunctionDescription(
    "SparseAccelerator", 
    "Dot",                            # SparseAccelerator.Dot(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const dot_Desc = FunctionDescription(
    "Main", 
    "dot",                            # SparseAccelerator.Dot(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const copy_Desc = FunctionDescription(
    "Main", 
    "copy",                           # 
    (Vector, ),                       # The arguments must be a vector
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1]))                   # The return value (0) and the argument are inter-dependent
)

const WAXPBY_Desc = FunctionDescription(
    "SparseAccelerator", 
    "WAXPBY",                         # SparseAccelerator.WAXPBY(alpha::Number, x::Vector, beta::Number, y::Vector)
    (Number, Vector, Number, Vector), 
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 2; 4]))                # The return value (0) and the two arguments (2, 4) are inter-dependent
)

const WAXPBY!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "WAXPBY!",                        # SparseAccelerator.WAXPBY!(w::Vector, alpha::Number, x::Vector, beta::Number, y::Vector)
    (Vector, Number, Vector, Number, Vector), 
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 3; 5]))             # The return value (0) and the arguments 1, 3, 5 are inter-dependent. In fact, 0 and 1 are the same array
)

const add_vector_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (VECTOR_OR_NUM, VECTOR_OR_NUM),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const add_matrix_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (MATRIX_OR_NUM, MATRIX_OR_NUM),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const sub_vector_Desc = FunctionDescription(
    "Main", 
    "-",                              
    (VECTOR_OR_NUM, VECTOR_OR_NUM),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments (1, 2) are inter-dependent
)

const norm_Desc = FunctionDescription(
    "Main", 
    "norm",                              
    (AbstractArray, ),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    IA_NONE                           # No inter-dependent arrays
)

const spones_Desc = FunctionDescription(
    "Main", 
    "spones",                              
    (AbstractSparseMatrix, ),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1]))                   # The return value (0) and the argument are inter-dependent
)

const size_Desc = FunctionDescription(
    "Main", 
    "size",                              
    (AbstractArray, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    IA_NONE                           # No inter-dependent arrays
)

const Array_Desc = FunctionDescription(
    "Main", 
    "Array",                          # Array(Any, Integer)
    (Any, Integer),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is actually not distributive, but it constructs an array whose values are unitialized, and thus we can treat it as distributive 
    IA_NONE                           # No inter-dependent arrays
)

const fill!_Desc = FunctionDescription(
    "Main", 
    "fill!",                          
    (AbstractArray, Number),
    Set(1),                           # argument 1 is updated
    true,                             # The function is distributive 
    IA_NONE                           # No inter-dependent arrays: the result array is filled with the same value, and thus independent of the original array's content
)

const max_Desc = FunctionDescription(
    "Main", 
    "max",                          
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    ia(Set([0, 1]))                   # The return value (0) and the argument are inter-dependent
)

const scale_Desc = FunctionDescription(
    "Main", 
    "scale",                          
    (AbstractMatrix, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments are inter-dependent
)

const sum_Desc = FunctionDescription(
    "Main", 
    "sum",                          
    (AbstractArray, ),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    IA_NONE                           # No inter-dependent arrays
)

const convert_Desc = FunctionDescription(
    "Main", 
    "convert",                          
    (Type{Array{Float64,1}}, Vector), # Cannot use Type{Vector} as Type{Array{Float64,1}} <: Type{Vector} is false! but Type{Array{Float64,1}} <: Type{Array{Float64,1}} is true! Not sure why
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    ia(Set([0; 2]))                   # The return value (0) and the argument 2 are inter-dependent
)

const eltype_Desc = FunctionDescription(
    "Main", 
    "eltype",                          
    (AbstractArray, ), 
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    IA_NONE                           # No inter-dependent arrays
)

const inverse_divide_Desc = FunctionDescription(
    "Main", 
    "\\",                              
    (AbstractSparseMatrix, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments are inter-dependent
)

const fwdTriSolve!_Desc = FunctionDescription(
    "Main", 
    "fwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments are inter-dependent
)

const bwdTriSolve!_Desc = FunctionDescription(
    "Main", 
    "bwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    ia(Set([0; 1; 2]))                # The return value (0) and the arguments are inter-dependent
)

function_descriptions  = [
    PointwiseMultiply_Desc,
    PointwiseMultiply!_Desc,
    PointwiseDivide_Desc,
    PointwiseDivide1_Desc,
    PointwiseDivide!_Desc,
    SpMV_Desc,
    star_Desc,
    star1_Desc,
    star2_Desc,
    Dot_Desc,
    dot_Desc,
    WAXPBY_Desc,
    WAXPBY!_Desc,
    add_vector_Desc,
    add_matrix_Desc,
    sub_vector_Desc,
    norm_Desc,
    spones_Desc,
    size_Desc,
    Array_Desc,
    fill!_Desc,
    max_Desc,
    scale_Desc,
    sum_Desc,
    convert_Desc,
    eltype_Desc,
    inverse_divide_Desc,
    copy_Desc,
    fwdTriSolve!_Desc,
    bwdTriSolve!_Desc
]

function show_function_descriptions()
    println("Function descriptions: ", function_descriptions)
end

function look_for_function_description(
    module_name    :: String, 
    function_name  :: String, 
    argument_types :: Tuple{Type}
)
    return look_for_function(function_descriptions, module_name, function_name, argument_types)
end
