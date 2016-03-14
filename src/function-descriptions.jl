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
const COLUMN_ROW         = 3 # A's column permutation vector is B's row permutation vector.

# A function is assumed to have the following form
#      f(argument1, argument2, ...)
# In the description, its return result is represented by number 0, 
# and its arguments by number 1, 2, and so on.
# Some arguments may be output as well, which will be described by
# the "output" field.
# If the return result is a tuple, the elements of the tuple are denoted 0.1,
# 0.2, etc.
# ASSUMPTION: all the described is "must" information. That is, the information
# is true in all cases.
# ISSUE: depending on the specific inputs, a function "may" have different
# behavior. For example, it may update some arguments, or not update them.
# Two arguments (arrays) may be inter-dependent in this case, or not in another
# case. The function may be distributive in this case, but not in another.
# How to represent and use some "may" information?

immutable FunctionDescription
    module_name     :: AbstractString # Module of the function. It is "nothing" if it is not in a Julia module, e.g. if this is a C function 
    function_name   :: AbstractString # Name of the function
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
          (0, 2, ROW_ROW),
          (1, 2, ROW_ROW)
    ])
)

const element_wise_multiply1_Desc = FunctionDescription(
    "Main", 
    ".*",
    (Vector, Vector),                 # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW),
          (1, 2, ROW_ROW)
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

@doc """ A*x """
const SpMV_2_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",                           # SparseAccelerator.SpMV(A::SparseMatrixCSC, x::Vector)
    (SparseMatrixCSC, Vector),        # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, COLUMN_ROW_INVERSE)
    ])
)

@doc """ alpha*A*x """
const SpMV_3_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",
    (Number, SparseMatrixCSC, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE)
    ])
)

@doc """ alpha*A*x + y """
const SpMV_4_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",
    (Number, SparseMatrixCSC, Vector, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE),
          (0, 4, ROW_ROW),
          (2, 4, ROW_ROW)
    ])
)

@doc """ alpha*A*x + beta*y """
const SpMV_5_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",
    (Number, SparseMatrixCSC, Vector, Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE),
          (0, 5, ROW_ROW),
          (2, 5, ROW_ROW)
    ])
)

@doc """ alpha*A*x + beta*y + gamma """
const SpMV_6_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV",
    (Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE),
          (0, 5, ROW_ROW),
          (2, 5, ROW_ROW)
    ])
)

@doc """ SparseAccelerator.SpMV!(A, x) """
const SpMV!_3_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV!",
    (Vector, SparseMatrixCSC, Vector),
    Set(1),                           # Argument 1 (the vector) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW),
          (1, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE)
    ])
)

@doc """ w = alpha*A*x + beta*y + gamma """
const SpMV!_7_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV!",
    (Vector, Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    Set(1),                           # Argument 1 (the vector) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 3, ROW_ROW),
          (0, 6, ROW_ROW),
          (1, 3, ROW_ROW),
          (1, 6, ROW_ROW),
          (3, 4, COLUMN_ROW_INVERSE),
          (3, 6, ROW_ROW)
    ])
)

@doc """ w = (alpha*A*x + beta*y + gamma) .* z """
const SpMV!_8_parameters_Desc = FunctionDescription(
    "SparseAccelerator", 
    "SpMV!",
    (Vector, Number, SparseMatrixCSC, Vector, Number, Vector, Number, Vector),
    Set(1),                           # Argument 1 (the vector) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 3, ROW_ROW),
          (0, 6, ROW_ROW),
          (0, 8, ROW_ROW),
          (1, 3, ROW_ROW),
          (1, 6, ROW_ROW),
          (1, 8, ROW_ROW),
          (3, 4, COLUMN_ROW_INVERSE),
          (3, 6, ROW_ROW),
          (3, 8, ROW_ROW),
          (6, 8, ROW_ROW)
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

const star3_Desc = FunctionDescription(
    "Main", 
    "*",
    (SparseMatrixCSC, SparseMatrixCSC),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, COLUMN_ROW_INVERSE),
          (0, 2, COLUMN_COLUMN)
    ])
)

const A_mul_B!_Desc = FunctionDescription(
    "Main", 
    "A_mul_B!",
    (Vector, SparseMatrixCSC, Vector),
    Set(1),                           # argument 1 (w) is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW),
          (1, 2, ROW_ROW),
          (2, 3, COLUMN_ROW_INVERSE)
    ])
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

const abs_Desc = FunctionDescription(
    "Main", 
    "abs",
    (Vector, ),                       # The arguments must be a vector
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

const log_Desc = FunctionDescription(
    "Main", 
    "log",
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

const exp_Desc = FunctionDescription(
    "Main", 
    "exp",                              
    (Vector, ),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
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

@doc """ max(A, value) """
const max_Desc = FunctionDescription(
    "Main", 
    "max",                          
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive 
    Set([ (0, 1, ROW_ROW) ])
)

@doc """ max(value, A) """
const max1_Desc = FunctionDescription(
    "Main", 
    "max",
    (Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW) ])
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
    Set([ (1, 2, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

const fwdTriSolve!1_Desc = FunctionDescription(
    "SparseAccelerator", 
    "fwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    Set([ (1, 2, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

const fwdTriSolve!2_Desc = FunctionDescription(
    "SparseAccelerator", 
    "fwdTriSolve!",                              
    (Vector, AbstractSparseMatrix, Vector),
    Set(1),                           # Argument 1 (the vector) is updated
    true,                             # The function is distributive
    Set([ (2, 1, COLUMN_ROW_INVERSE),
          (2, 3, ROW_ROW)
    ])
)

const bwdTriSolve!_Desc = FunctionDescription(
    "Base.SparseMatrix", 
    "bwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    Set([ (1, 2, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

const bwdTriSolve!1_Desc = FunctionDescription(
    "SparseAccelerator", 
    "bwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),                           # Argument 2 (the vector) is updated
    true,                             # The function is distributive
    Set([ (1, 2, COLUMN_ROW_INVERSE),
          (1, 2, ROW_ROW)
    ])
)

const bwdTriSolve!2_Desc = FunctionDescription(
    "SparseAccelerator", 
    "bwdTriSolve!",                              
    (Vector, AbstractSparseMatrix, Vector),
    Set(1),                           # Argument 1 (the vector) is updated
    true,                             # The function is distributive
    Set([ (2, 1, COLUMN_ROW_INVERSE),
          (2, 3, ROW_ROW)
    ])
)

const negative_Desc = FunctionDescription(
    "Main", 
    "-",
    (Vector,),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

const divide_Desc = FunctionDescription(
    "Main", 
    "/",
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

@doc """ A*x """
const Ac_mul_B_Desc = FunctionDescription(
    "Main", 
    "Ac_mul_B",
    (SparseMatrixCSC, Vector),        # The arguments must be vectors
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (1, 2, COLUMN_ROW_INVERSE)
    ])
)

@doc """ min(value, A) """
const min_Desc = FunctionDescription(
    "Main", 
    "min",
    (Number, Vector),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 2, ROW_ROW) ])
)

@doc """ min(A, value) """
const min1_Desc = FunctionDescription(
    "Main", 
    "min",
    (Vector, Number),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW) ])
)

#The following are made for context-test3 with reordering on. There were
#call resolution issues. For setfield!, there is also a question how to
#handle fields (field-sensitive analysis?).
 
#@doc """ typeof """
#const typeof_Desc = FunctionDescription(
#    "", 
#    "typeof",
#    (Any, ),
#    UPDATED_NONE,                     # No argument is updated
#    true,                             # The function is distributive
#    Set()
#)

#@doc """ convert """
#const convert_Desc = FunctionDescription(
#    "", 
#    "convert",
#    (Type{Array{Float64,1}}, Array{Float64,1}),
#    UPDATED_NONE,                     # No argument is updated
#    true,                             # The function is distributive
#    Set([ (0, 2, ROW_ROW) ])
#)

#@doc """ fieldtype """
#const fieldtype_Desc = FunctionDescription(
#    "", 
#    "fieldtype",
#    (Type{SparseMatrixCSC{Float64,Int32}}, QuoteNode),
#    UPDATED_NONE,                     # No argument is updated
#    true,                             # The function is distributive
#    Set()
#)

#@doc """ setfield!(A, nzval) = x """
#const setfield!_Desc = FunctionDescription(
#    "", 
#    "setfield!",
#    (Type{SparseMatrixCSC{Float64,Int32}}, QuoteNode(:nzval), Array{Float64,1}),
#    Set(1),
#    true,                             # The function is distributive
#    Set([ (1, 3, ROW_ROW) ])
#)

# Make assignment a special function
const asignment_Desc = FunctionDescription(
    "", 
    ":=",
    (Any, AbstractMatrix),            # If LHS is a GenSym, lamda may or may not have its type info (in Julia 0.4 rc1). So use Any.
    Set(1),                           # Argument 1 (the left hand side) is updated
    true,                             # The function is distributive
    Set([ (1, 2, ROW_ROW),
          (1, 2, COLUMN_COLUMN)
    ])
)

const asignment1_Desc = FunctionDescription(
    "", 
    ":=",
    (Any, Vector),                    # If LHS is a GenSym, lamda may or may not have its type info (in Julia 0.4 rc1). So use Any.
    Set(1),                           # Argument 1 (the left hand side) is updated
    true,                             # The function is distributive
    Set([ (1, 2, ROW_ROW) ])
)

# Make apply_type as a special assignment
const apply_type_Desc = FunctionDescription(
    "", 
    "apply_type",
    (AbstractMatrix,),                # If LHS is a GenSym, lamda may or may not have its type info (in Julia 0.4 rc1). So use Any.
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN)
    ])
)

const Printf_Desc = FunctionDescription(
    "Base", 
    "Printf",
    (Any, Any),                    # If LHS is a GenSym, lamda may or may not have its type info (in Julia 0.4 rc1). So use Any.
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    IA_NONE                           # No inter-dependent arrays
)

# The following 3 are specific to LBFGS
const lbfgs_compute_direction!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "lbfgs_compute_direction!",
    (Array{Float64,1}, Int64, Int64, Int64, Array{Float64,2}, Array{Float64,2}, Array{Float64,1}),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (1, 5, ROW_ROW),
          (1, 6, ROW_ROW),
          (1, 7, ROW_ROW)
    ])
)

const lbfgs_loss_function1_Desc = FunctionDescription(
    "SparseAccelerator", 
    "lbfgs_loss_function1",
    (Vector, Vector, Number),         # If LHS is a GenSym, lamda may or may not have its type info (in Julia 0.4 rc1). So use Any.
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    IA_NONE                           # No inter-dependent arrays
)


const lbfgs_loss_function2!_Desc = FunctionDescription(
    "SparseAccelerator", 
    "lbfgs_loss_function2!",
    (Vector, Vector, Vector),         # If LHS is a GenSym, lamda may or may not have its type info (in Julia 0.4 rc1). So use Any.
    Set(1),
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 2, ROW_ROW),
          (1, 2, ROW_ROW)
    ])
)

const ilu_Desc = FunctionDescription(
    "SparseAccelerator", 
    "ilu",
    (SparseMatrixCSC,),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ ("0.1", 1, ROW_ROW),          # The return value is a tuple. Both elements
          ("0.1", 1, COLUMN_COLUMN),    #   in the tuple are inter-dependent with the
          ("0.2", 1, ROW_ROW),          #   input.
          ("0.2", 1, COLUMN_COLUMN),
    ])
)

const tril_Desc = FunctionDescription(
    "Main", 
    "tril",
    (SparseMatrixCSC,),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN),
    ])
)

const triu_Desc = FunctionDescription(
    "Main", 
    "triu",
    (SparseMatrixCSC,),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN),
    ])
)

const diag_Desc = FunctionDescription(
    "Main", 
    "diag",
    (SparseMatrixCSC,),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),
          (0, 1, COLUMN_COLUMN),
    ])
)

const spdiagm_Desc = FunctionDescription(
    "Main", 
    "spdiagm",
    (Vector,),
    UPDATED_NONE,                     # No argument is updated
    true,                             # The function is distributive
    Set([ (0, 1, ROW_ROW),            # Both the rows and the columns of the result
          (0, 1, COLUMN_ROW),         # are dependent on the vector
    ])
)

function_descriptions  = [
    element_wise_multiply_Desc,
    element_wise_multiply1_Desc,
    element_wise_multiply!_Desc,
    element_wise_divide_Desc,
    element_wise_divide1_Desc,
    element_wise_divide2_Desc,
    element_wise_divide!_Desc,
    SpMV_2_parameters_Desc,
    SpMV_3_parameters_Desc,
    SpMV_4_parameters_Desc,
    SpMV_5_parameters_Desc,
    SpMV_6_parameters_Desc,
    SpMV!_3_parameters_Desc,
    SpMV!_7_parameters_Desc,
    SpMV!_8_parameters_Desc,
    star_Desc,
    star1_Desc,
    star2_Desc,
    star3_Desc,
    A_mul_B!_Desc,
    Dot_Desc,
    dot_Desc,
    copy_Desc,
    abs_Desc,
    log_Desc,
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
    exp_Desc,
    norm_Desc,
    spones_Desc,
    size_Desc,
    Array_Desc,
    fill!_Desc,
    max_Desc,
    max1_Desc,
    scale_Desc,
    sum_Desc,
    convert_Desc,
    eltype_Desc,
    inverse_divide_Desc,
    fwdTriSolve!_Desc,
    fwdTriSolve!1_Desc,
    fwdTriSolve!2_Desc,
    bwdTriSolve!_Desc,
    bwdTriSolve!1_Desc,
    bwdTriSolve!2_Desc,
    negative_Desc,
    divide_Desc,
    Ac_mul_B_Desc,
    min_Desc,
    min1_Desc,
    asignment_Desc,
    asignment1_Desc,
    apply_type_Desc,
    lbfgs_compute_direction!_Desc,
    lbfgs_loss_function1_Desc,
    lbfgs_loss_function2!_Desc,
    ilu_Desc,
    tril_Desc,
    triu_Desc,
    diag_Desc,
    spdiagm_Desc
]

function look_for_function_description(
    module_name    :: AbstractString, 
    function_name  :: AbstractString, 
    argument_types :: Tuple
)
    return look_for_function(function_descriptions, module_name, function_name, argument_types)
end
