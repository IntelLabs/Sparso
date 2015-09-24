# This file contains a description for a function, including the function's
# input/output arguments, distributivity, and inter-dependent arrays. 
# This gives compiler accurate information about the function without analyzing
# it. User can easily describe his/her own functions by adding entries to the
# vector below.

# Below are the possible relationships between two arrays in term of reordering.
# Assume the two arrays are A and B. Each matrix can have two permutation vectors
# with it, one for rows, one for columns, to reorder the matrix; and each
# vector can have one row permutation vector, to reorder the vector.
const ROW_ROW            = 0 # A's row permutation vector is the same as B's row permutation vector.
const COLUMN_COLUMN      = 1 # A's column permutation vector is the same as B's column permutation vector.
const COLUMN_ROW_INVERSE = 2 # A's column permutation vector is the inverse of B's row permutation vector.

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
    IA              :: Set    # Inter-dependence arrays. Each element is a triple
                              # (a, b, c), where a and b are two arguments' indices,
                              # and c is their relationship in terms of reordering.
                              # E.g. Set((0, 1, ROW_ROW), (3, 4, COLUMN_ROW_INVERSE)])
                              # means that result (0) and argument 1 are 
                              # inter-dependent, and their row permutation 
                              # vectors are the same; argument 3 and 4 are
                              # inter-dependent, and 3's column permutation vector
                              # is the inverse of the 4's row permutation vector.
end

const UPDATED_NONE  = Set()
const IA_NONE       = Set()

const element_wise_multiply_Desc = FunctionDescription(
    "SparseAccelerator", 
    "element_wise_multiply",          # SparseAccelerator.element_wise_multiply(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW)
         # No need to describe 1 and 2's relationship: the relationship is transitive.
    ])
)

const element_wise_multiply!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "element_wise_multiply!",         # SparseAccelerator.element_wise_multiply!(w::Vector, x::Vector, y::Vector)
    (Vector, Vector, Vector),         # The arguments must be vectors
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, ROW_ROW),
          (1, 3, ROW_ROW)
    ])
)

const element_wise_divide_Desc = FunctionDescription(
    "SparseAccelerator", 
    "element_wise_divide",            # SparseAccelerator.element_wise_divide(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW)
    ])
)

const element_wise_divide1_Desc = FunctionDescription(
    "Main", 
    "./",           
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW)
    ])
)

const element_wise_divide2_Desc = FunctionDescription(
    "Main", 
    "./",           
    (Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW) ])
)

const element_wise_divide!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "element_wise_divide!",           # SparseAccelerator.element_wise_divide!(w::Vector, x::Vector, y::Vector)
    (Vector, Vector, Vector),         # The arguments must be vectors
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, ROW_ROW),
          (1, 3, ROW_ROW)
    ])
)

const SpMV_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",                           # SparseAccelerator.SpMV(A::SparseMatrixCSC, x::Vector)
    (SparseMatrixCSC, Vector),        # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, COLUMN_ROW_INVERSE)
    ])
)

const star_Desc = FunctionDescription(
    "Main", 
    "*",                              # *(A::SparseMatrixCSC, x::Vector)
    (SparseMatrixCSC, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, COLUMN_ROW_INVERSE)
    ])
)

const star1_Desc = FunctionDescription(
    "Main", 
    "*",                              
    (Number, SparseMatrixCSC, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE)
    ])
)

const star2_Desc = FunctionDescription(
    "Main", 
    "*",           
    (Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW) ])
)

const Dot_Desc = FunctionDescription(
    "SparseAccelerator", 
    "dot",                            # SparseAccelerator.Dot(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (1, 2, ROW_ROW) ])
)

const dot_Desc = FunctionDescription(
    "Main", 
    "dot",                            # SparseAccelerator.Dot(x::Vector, y::Vector)
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (1, 2, ROW_ROW) ])
)

const copy_Desc = FunctionDescription(
    "Main", 
    "copy",
    (Vector, ),                       # The arguments must be a vector
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

const WAXPBY_Desc = FunctionDescription(
    "SparseAccelerator", 
    "WAXPBY",                         # SparseAccelerator.WAXPBY(alpha::Number, x::Vector, beta::Number, y::Vector)
    (Number, Vector, Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (2, 4, ROW_ROW)
    ])
)

const WAXPBY!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "WAXPBY!",                        # SparseAccelerator.WAXPBY!(w::Vector, alpha::Number, x::Vector, beta::Number, y::Vector)
    (Vector, Number, Vector, Number, Vector),
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 3, ROW_ROW),
          (3, 5, ROW_ROW)
    ])
)

const add_vector_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (Vector, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, ROW_ROW)
    ])
)

const add_vector1_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW) ])
)

const add_vector2_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

const sub_vector_Desc = FunctionDescription(
    "Main", 
    "-",                              
    (Vector, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, ROW_ROW)
    ])
)

const sub_vector1_Desc = FunctionDescription(
    "Main", 
    "-",                              
    (Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW) ])
)

const sub_vector2_Desc = FunctionDescription(
    "Main", 
    "-",                              
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

const add_matrix_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (AbstractMatrix, AbstractMatrix),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN),
          (1, 2, ROW_ROW),
          (1, 2, COLUMN_COLUMN)
    ])
)

const add_matrix1_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (Number, AbstractMatrix),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (0, 2, COLUMN_COLUMN)
    ])
)

const add_matrix2_Desc = FunctionDescription(
    "Main", 
    "+",                              
    (AbstractMatrix, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN)
    ])
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
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN)
    ])
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
    (AbstractMatrix, Number),
    Set(1),                           # argument 1 is updated
    true,                             # The function is distributive 
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN)
    ])
)

const max_Desc = FunctionDescription(
    "Main", 
    "max",                          
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    Set([ (0, 1, ROW_ROW) ])
)

const scale_Desc = FunctionDescription(
    "Main", 
    "scale",                          
    (AbstractMatrix, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    Set([ (0, 1, ROW_ROW),
          (1, 2, COLUMN_ROW_INVERSE)
    ])
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
    (Array{Float64,1}, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    Set([ (0, 2, ROW_ROW) ])
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
    Set([ (1, 0, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

const fwdTriSolve!_Desc = FunctionDescription(
    "Base.SparseMatrix", 
    "fwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    Set([ (1, 0, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

const bwdTriSolve!_Desc = FunctionDescription(
    "Base.SparseMatrix", 
    "bwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    Set([ (1, 0, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

# Make assignment a special function
const asignment_Desc = FunctionDescription(
    "", 
    ":=",
    (AbstractSparseMatrix, AbstractSparseMatrix),
    Set(1),                           # Argument 1 (the left hand side) is updated
    true,                             # The function is distributive
    Set([ (1, 2, ROW_ROW),
          (1, 2, COLUMN_COLUMN)
    ])
)

const asignment1_Desc = FunctionDescription(
    "", 
    ":=",
    (Vector, Vector),
    Set(1),                           # Argument 1 (the left hand side) is updated
    true,                             # The function is distributive
    Set([ (1, 2, ROW_ROW) ])
)

function_descriptions  = [
    element_wise_multiply_Desc,
    element_wise_multiply!_Desc,
    element_wise_divide_Desc,
    element_wise_divide1_Desc,
    element_wise_divide2_Desc,
    element_wise_divide!_Desc,
    SpMV_Desc,
    star_Desc,
    star1_Desc,
    star2_Desc,
    Dot_Desc,
    dot_Desc,
    copy_Desc,
    WAXPBY_Desc,
    WAXPBY!_Desc,
    add_vector_Desc,
    add_vector1_Desc,
    add_vector2_Desc,
    sub_vector_Desc,
    sub_vector1_Desc,
    sub_vector2_Desc,
    add_matrix_Desc,
    add_matrix1_Desc,
    add_matrix2_Desc,
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
    fwdTriSolve!_Desc,
    bwdTriSolve!_Desc,
    asignment_Desc,
    asignment1_Desc
]

function look_for_function_description(
    module_name    :: String, 
    function_name  :: String, 
    argument_types :: Tuple
)
    return look_for_function(function_descriptions, module_name, function_name, argument_types)
end