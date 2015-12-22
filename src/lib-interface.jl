# This file contains the routines that interfaces with SpMP library.

const LIB_PATH = libcsr

@doc """
Create a knob for a matrix with the given properties of it, while this matrix
might not have been created yet.

Properties: 
constant_valued:        The matrix is a constant in value. It implies 
                        constant_structured below.
constant_structured:    The matrix has always the same structure, even if its
                        value may be changed.
is_symmetric:           The matrix is symmetric in value.
                        It implies is_structure_symmetric below.
is_structure_symmetric: The matrix is symmetric in structure.
is_structure_only :     Only the structure of matrix should be used.
is_single_def:          The matrix is statically defined only once.
"""
function new_matrix_knob(
    A                      :: Symbol, # Useless. Only for better readability of the generated code
    constant_valued        = false,
    constant_structured    = false,
    is_symmetric           = false,
    is_structure_symmetric = false,
    is_structure_only      = false,
    is_single_def          = false
 )
    # Note: is_structure_symmetric is only for CoSP2. In general, the matrix
    # must be constant valued or structured.
    assert(constant_valued || constant_structured || is_structure_symmetric || !collective_structure_prediction_enabled)
    assert(!constant_valued || constant_structured || !collective_structure_prediction_enabled)
    assert(!is_symmetric || is_structure_symmetric || !collective_structure_prediction_enabled)
    mknob = ccall((:NewMatrixKnob, LIB_PATH), Ptr{Void},
                   (Bool, Bool, Bool, Bool, Bool, Bool),
                   constant_valued, constant_structured, is_symmetric,
                   is_structure_symmetric, is_structure_only, is_single_def)
end

@doc """ Fill matrix info into the knob """
function fill_matrix_info_once(
    mknob :: Ptr{Void},
    A     :: SparseMatrixCSC
 )
    ccall((:FillMatrixInfoOnce, LIB_PATH), Void,
            (Ptr{Void}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
             mknob, A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval))
end

@doc """ Set the matrix as constant_structured. """
function set_constant_valued(
    mknob :: Ptr{Void}
)
    ccall((:SetConstantValued, LIB_PATH), Void, (Ptr{Void},), mknob)
end

@doc """ Set the matrix as constant_structured. """
function set_constant_structured(
    mknob :: Ptr{Void}
)
    ccall((:SetConstantStructured, LIB_PATH), Void, (Ptr{Void},), mknob)
end

@doc """ Set the matrix as value symmetric. """
function set_value_symmetric(
    mknob :: Ptr{Void}
)
    ccall((:SetValueSymmetric, LIB_PATH), Void, (Ptr{Void},), mknob)
end

@doc """ Set the matrix as structurally symmetric. """
function set_structure_symmetric(
    mknob :: Ptr{Void}
)
    ccall((:SetStructureSymmetric, LIB_PATH), Void, (Ptr{Void},), mknob)
end

@doc """ set matrix derivative """
function set_derivative(
    mknob           :: Ptr{Void},
    derivative_type :: Int,
    derivative      :: Ptr{Void}
)
    ccall((:SetDerivative, LIB_PATH), Void, (Ptr{Void}, Cint, Ptr{Void}), mknob, derivative_type, derivative)
end

@doc """
Derivative types between two matrix knobs.
NOTE: when adding any new type, or changing any type' integer value,
make sure to change int2derivative_map as well.
"""
const DERIVATIVE_TYPE_TRANSPOSE = 0
const DERIVATIVE_TYPE_SYMMETRIC = 1
const DERIVATIVE_TYPE_LOWER_TRIANGULAR = 2
const DERIVATIVE_TYPE_UPPER_TRIANGULAR = 3

const int2derivative_map = Dict(
    0 => GlobalRef(SparseAccelerator, :DERIVATIVE_TYPE_TRANSPOSE),
    1 => GlobalRef(SparseAccelerator, :DERIVATIVE_TYPE_SYMMETRIC),
    2 => GlobalRef(SparseAccelerator, :DERIVATIVE_TYPE_LOWER_TRIANGULAR),
    3 => GlobalRef(SparseAccelerator, :DERIVATIVE_TYPE_UPPER_TRIANGULAR)
)

@doc """ Delete a matrix knob. """
function delete_matrix_knob(
    mknob :: Ptr{Void}
)
    ccall((:DeleteMatrixKnob, LIB_PATH), Void, (Ptr{Void},), mknob)
end

@doc """ Propagate a matrix knob's information to another. """
function propagate_matrix_info(
    to_mknob   :: Ptr{Void},
    from_mknob :: Ptr{Void}
)
    ccall((:PropagateMatrixInfo, LIB_PATH), Void, (Ptr{Void}, Ptr{Void}), to_mknob, from_mknob)
end

@doc """ Associate a matrix knob with a function knob. """
function add_mknob_to_fknob(
    mknob :: Ptr{Void}, 
    fknob :: Ptr{Void}
)
    ccall((:AddMatrixKnob, LIB_PATH), Void, (Ptr{Void}, Ptr{Void}), 
        fknob, mknob)
end

@doc """ Create a function knob. """
function new_function_knob(
    fknob_creator :: AbstractString
)
    assert(fknob_creator != "")

    # Simulate ccall((fknob_creator, LIB_PATH), Ptr{Void}, ()). We cannot directly
    # use this ccall here because (fknob_creator, LIB_PATH) is treated as a tuple,
    # instead of a pointer or expression.
    expr = Expr(:call, 
        TopNode(:ccall), 
        Expr(:call, TopNode(:tuple), QuoteNode(fknob_creator), LIB_PATH), 
        Ptr{Void}, 
        Expr(:call, TopNode(:svec))
    )
    eval(expr)
end
new_function_knob() = new_function_knob("NewFunctionKnob")

@doc """ Delete a function knob. """
function delete_function_knob(
    fknob_deletor :: AbstractString, 
    fknob         :: Ptr{Void}
)
    assert(fknob_deletor != "")

    # Simulate ccall((fknob_deletor, LIB_PATH), Void, (Ptr{Void},), fknob). We cannot directly
    # use this ccall here because (fknob_deletor, LIB_PATH) is treated as a tuple,
    # instead of a pointer or expression.
    expr = Expr(:call, 
        TopNode(:ccall), 
        Expr(:call, TopNode(:tuple), QuoteNode(fknob_deletor), LIB_PATH), 
        GlobalRef(Main, :Void), 
        Expr(:call, 
                TopNode(:svec),
                Expr(:call, 
                      TopNode(:apply_type), 
                      GlobalRef(Main, :Ptr),
                      GlobalRef(Main, :Void)
                )
        ),
        Expr(:call, 
                TopNode(:unsafe_convert), 
                Expr(:call, 
                      TopNode(:apply_type), 
                      GlobalRef(Main, :Ptr),
                      GlobalRef(Main, :Void)
                ),
                fknob
        ),
        fknob
    )
    eval(expr)
end
delete_function_knob(fknob :: Ptr{Void}) = delete_function_knob("DeleteFunctionKnob", fknob)

function replace_lower_with_UT(
    A     :: SparseMatrixCSC,
    U     :: SparseMatrixCSC
)
    ccall((:ReplaceLowerWithUT, LIB_PATH), Void,
           (
            Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}
           ),
            A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
            U.m, U.n, pointer(U.colptr), pointer(U.rowval), pointer(U.nzval))
end

@doc """
Context-sensitive forward triangular solver, equivalent to Base.SparseMatrix.
fwdTriSolve!(L, b).
"""
function fwdTriSolve!(
    L     :: SparseMatrixCSC{Float64, Int32}, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    mknob_L = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    assert(mknob_L != C_NULL)
    fill_matrix_info_once(mknob_L, L)

    ccall((:ForwardTriangularSolve, LIB_PATH), Void,
           (
            Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void},
           ),
            L.m, L.n, pointer(L.colptr), pointer(L.rowval), pointer(L.nzval),
            pointer(b), pointer(b), fknob)
    b
end

@doc """
Context-sensitive forward triangular solver, equivalent to x = Base.SparseMatrix.
fwdTriSolve!(L, b).
"""
function fwdTriSolve!(
    x     :: Vector,
    L     :: SparseMatrixCSC{Float64, Int32}, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    # ISSUE: if x is not allocated space yet in the caller,can we really copy b to x in this way?
    x[:] = copy(b)
    fwdTriSolve!(L, x, fknob)
    x
end

@doc """ 
Context-sensitive backward triangular solver. 
"""
function bwdTriSolve!(
    U     :: SparseMatrixCSC{Float64, Int32}, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    mknob_U = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    assert(mknob_U != C_NULL)
    fill_matrix_info_once(mknob_U, U)
    
    ccall((:BackwardTriangularSolve, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
               U.m, U.n, pointer(U.colptr), pointer(U.rowval), pointer(U.nzval),
               pointer(b), pointer(b), fknob)
    b
end

@doc """
Context-sensitive backward triangular solver, equivalent to x = Base.SparseMatrix.
bwdTriSolve!(L, b).
"""
function bwdTriSolve!(
    x     :: Vector,
    U     :: SparseMatrixCSC{Float64, Int32}, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    # ISSUE: if x is not allocated space yet in the caller,can we really copy b to x in this way?
    x[:] = copy(b)
    bwdTriSolve!(U, x, fknob)
    x
end

function speye_int32(m::Integer)
  rowval = [Int32(x) for x in [1:m;]]
  colptr = [Int32(x) for x in [rowval; m + 1]]
  nzval = ones(Float64, m)
  return SparseMatrixCSC(m, m, colptr, rowval, nzval)
end

function ilu(
    A     :: SparseMatrixCSC{Float64, Int32}
 )
    LU_nzval = zeros(nnz(A))

    ccall((:ILU, LIB_PATH), Void,
          (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
          size(A, 1), size(A, 2),
          A.colptr, A.rowval, A.nzval,
          LU_nzval,
          C_NULL)

    LU = SparseMatrixCSC{Cdouble, Cint}(size(A, 1), size(A, 2), copy(A.colptr), copy(A.rowval), LU_nzval)
    LU = LU'

    L = tril(LU, -1) + speye_int32(size(LU, 1))
    U = triu(LU)

    L, U
end

@doc """ 
Reorder a sparse matrix A and store the result in new_A. A itself is not changed.
If get_permutation is true, then compute the permutation and inverse permutation
vector P and inverse_P. Otherwise, P and inverse_P are given inputs.
One_based_input/output tells if A and new_A are 1-based.
ISSUE: We hard code a sparse matrix format to be SparseMatrixCSC{Cdouble, Cint},
and a permutation and inverse permutation vector to be Vector{Cint}.
Should make it general in future
"""
function reorder_matrix(
    A                :: SparseMatrixCSC{Float64, Int32}, 
    new_A            :: SparseMatrixCSC{Float64, Int32}, 
    P                :: Vector, 
    inverse_P        :: Vector
)
    ccall((:ReorderMatrix, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(new_A.colptr), pointer(new_A.rowval), pointer(new_A.nzval),
               pointer(P), pointer(inverse_P))
end

function reorder_matrix!(
    A                :: SparseMatrixCSC{Float64, Int32}, 
    P                :: Vector, 
    inverse_P        :: Vector,
)
    ccall((:ReorderMatrixInplace, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(P), pointer(inverse_P))
    A
end

@doc """ Inversely reorder vector V with permutation vector P """
function inverse_reorder_vector!(
    V     :: Vector, 
    P     :: Vector
)
    ccall((:reorderVectorWithInversePermInplace, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(P), length(V))
    V
end

function set_reordering_decision_maker(
    fknob :: Ptr{Void}
 )
    ccall((:SetReorderingDecisionMaker, LIB_PATH), Void,
         (Ptr{Void},),
         fknob)
end

function get_reordering_vectors(
    fknob :: Ptr{Void}
 )
    m = Ref{Cint}(0)
    n = Ref{Cint}(0)
    row_perm = ccall((:GetRowReorderingVector, LIB_PATH), Ptr{Cint},
         (Ptr{Void}, Ref{Cint}),
         fknob, m)
    row_inv_perm = ccall((:GetRowInverseReorderingVector, LIB_PATH), Ptr{Cint},
         (Ptr{Void}, Ref{Cint}),
         fknob, m)
    col_perm = ccall((:GetColReorderingVector, LIB_PATH), Ptr{Cint},
         (Ptr{Void}, Ref{Cint}),
         fknob, n)
    col_inv_perm = ccall((:GetColInverseReorderingVector, LIB_PATH), Ptr{Cint},
         (Ptr{Void}, Ref{Cint}),
         fknob, n)
    pointer_to_array(row_perm, m[]), pointer_to_array(row_inv_perm, m[]), pointer_to_array(col_perm, n[]), pointer_to_array(col_inv_perm, n[])
end

@doc """
Return a new representation of the sparse matrix, which is in CSC 
format, in CSR format.
"""
function create_CSR(A :: SparseMatrixCSC)
    ccall((:CSR_Create, LIB_PATH), Ptr{Void},
        (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        A.n, A.m, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval))
end

@doc """ Destroy the CSR representation of the sparse matrix. """
function destroy_CSR(A :: Ptr{Void})
    ccall((:CSR_Destroy, LIB_PATH), Void, (Ptr{Void},), A)
end

@doc """ w = (alpha*A*x + beta*y + gamma).*z """
function SpMV!(
    w     :: Vector, 
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector, 
    gamma :: Number,
    z     :: Vector,
    fknob :: Ptr{Void} = C_NULL # fknob == C_NULL => A is symmetric so we can directly call CSR_MultiplyWithVector with CSC matrix
)
    assert(length(w) == length(y) || length(y) == 0)
    assert(length(w) == length(z) || length(z) == 0)

    if use_SPMP
        if fknob == C_NULL
          A1 = create_CSR(A)
          ccall((:CSR_MultiplyWithVector, LIB_PATH), Void,
                (Ptr{Cdouble}, Cdouble, Ptr{Void}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
                pointer(w), alpha, A1, pointer(x), beta, length(y) == 0 ? C_NULL : pointer(y), gamma, length(z) == 0 ? C_NULL : pointer(z))
          destroy_CSR(A1)
        else
          mknob_A = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
          assert(mknob_A != C_NULL)
          fill_matrix_info_once(mknob_A, A)

          ccall((:SpMV, LIB_PATH), Void,
                (Cint, Cint,
                 Ptr{Cdouble},
                 Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
                 Ptr{Cdouble},
                 Cdouble, Ptr{Cdouble},
                 Cdouble,
                 Ptr{Cdouble},
                 Ptr{Void}),
                size(A, 1), size(A, 2),
                w,
                alpha, A.colptr, A.rowval, A.nzval,
                x,
                beta, length(y) == 0 ? C_NULL : y,
                gamma,
                length(z) == 0 ? C_NULL : z,
                fknob)
        end
        w
    else
        # use Julia implementation
        w[:] = alpha * A * x + beta * y + gamma
    end
end

@doc """ w = alpha*A*x + beta*y + gamma """
SpMV!(
    w     :: Vector, 
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector, 
    gamma :: Number,
    fknob :: Ptr{Void} = C_NULL) =
    SpMV!(w, alpha, A, x, beta, y, gamma, Array{Float64,1}[], fknob)

@doc """ y = A*x """
SpMV!(y :: Vector, A :: SparseMatrixCSC, x :: Vector, fknob :: Ptr{Void} = C_NULL) = 
    SpMV!(y, one(eltype(x)), A, x, zero(eltype(y)), y, zero(eltype(x)), fknob)

@doc """ w = alpha*A*x + beta*y """
function SpMV!(
    w     :: Vector{Complex128}, 
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector{Complex128}, 
    beta  :: Number, 
    y     :: Vector{Complex128}, 
    fknob :: Ptr{Void} = C_NULL # fknob == C_NULL => A is symmetric so we can directly call CSR_MultiplyWithVector with CSC matrix
)
    assert(length(w) == length(y) || length(y) == 0)

    if use_SPMP
        if fknob == C_NULL
          A1 = create_CSR(A)
          ccall((:CSR_MultiplyWithComplexVector, LIB_PATH), Void,
                (Ptr{Complex128}, Cdouble, Ptr{Void}, Ptr{Complex128}, Cdouble, Ptr{Complex128}, Cdouble, Ptr{Complex128}),
                pointer(w), alpha, A1, pointer(x), beta, length(y) == 0 ? C_NULL : pointer(y), 0, C_NULL)
          destroy_CSR(A1)
        else
          mknob_A = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
          assert(mknob_A != C_NULL)
          fill_matrix_info_once(mknob_A, A)

          ccall((:SpMVComplex, LIB_PATH), Void,
                (Cint, Cint,
                 Ptr{Complex128},
                 Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
                 Ptr{Complex128},
                 Cdouble, Ptr{Complex128},
                 Cdouble,
                 Ptr{Complex128},
                 Ptr{Void}),
                size(A, 1), size(A, 2),
                w,
                alpha, A.colptr, A.rowval, A.nzval,
                x,
                beta, length(y) == 0 ? C_NULL : y,
                0,
                C_NULL,
                fknob)
        end
    else
        # use Julia implementation
        w[:] = alpha * A * x + beta * y
    end
end

@doc """ w = alpha*A*x + beta*y + gamma """
# TODO: need a real 8 args SpMV! gamma is not handled yet.
SpMV!(
    w     :: Vector{Complex128}, 
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector{Complex128}, 
    beta  :: Number, 
    y     :: Vector{Complex128}, 
    gamma :: Number,
    fknob :: Ptr{Void} = C_NULL) =
    SpMV!(w, alpha, A, x, beta, y, fknob)

@doc """ alpha*A*x + beta*y + gamma """
function SpMV(
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector, 
    gamma :: Number,
    fknob :: Ptr{Void} = C_NULL
)
   # When fknob is not null, we just use
   # the CSC input matrix as if it's CSR
   # assuming the input matrix is symmetric
   w = Array(Cdouble, size(A, 1))
   SpMV!(w, alpha, A, x, beta, y, gamma, fknob)
   w
end

@doc """ alpha*A*x + beta*y """
SpMV(alpha :: Number, A :: SparseMatrixCSC, x :: Vector, beta :: Number, y :: Vector, fknob :: Ptr{Void} = C_NULL) =
    SpMV(alpha, A, x, beta, y, zero(eltype(x)), fknob)

@doc """ alpha*A*x + y """
SpMV(alpha :: Number, A :: SparseMatrixCSC, x :: Vector, y :: Vector, fknob :: Ptr{Void} = C_NULL) = 
    SpMV(alpha, A, x, one(eltype(y)), y, zero(eltype(x)), fknob)

@doc """ alpha*A*x """
SpMV(alpha :: Number, A :: SparseMatrixCSC, x :: Vector, fknob :: Ptr{Void} = C_NULL) = 
    SpMV(alpha, A, x, zero(eltype(x)), Array{Float64,1}[], zero(eltype(x)), fknob)

@doc """ A*x """
SpMV(A :: SparseMatrixCSC, x :: Vector, fknob :: Ptr{Void} = C_NULL) = 
    SpMV(one(eltype(x)), A, x, zero(eltype(x)), Array{Float64,1}[], zero(eltype(x)), fknob)

@doc """ Dot product of vector x and y """
function dot(
    x :: Vector, 
    y :: Vector
)
  assert(length(x) == length(y))

  if use_SPMP #&& length(x) > 4096
    ccall((:dot, LIB_PATH), Cdouble,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(x), pointer(y))
  else
    Base.dot(x, y)
  end
end

@doc """ Dot product of vector x and y """
function dot(
    x :: Vector{Complex128}, 
    y :: Vector{Complex128}
)
  assert(length(x) == length(y))

  if use_SPMP #&& length(x) >= 3072
    ccall((:dot, LIB_PATH), Complex128,
          (Cint, Ptr{Complex128}, Ptr{Complex128}),
          length(x), pointer(x), pointer(y))
  else
    Base.dot(x, y)
  end
end

@doc """ Dot product of vector x and y """
function dot(
    x :: Vector{Complex128}, 
    y :: Vector{Complex128}
)
  assert(length(x) == length(y))

  if use_SPMP
    ccall((:dot, LIB_PATH), Complex128,
          (Cint, Ptr{Complex128}, Ptr{Complex128}),
          length(x), pointer(x), pointer(y))
  else
    dot(x, y)
  end
end

@doc """ 2-norm of vector x """
function norm(
    x :: Vector
)
  sqrt(dot(x, x))
end

@doc """ 2-norm of vector x """
function norm(
    x :: Vector{Complex128}
)

  if use_SPMP #&& length(x) >= 8192
    ccall((:norm_complex, LIB_PATH), Cdouble,
          (Cint, Ptr{Complex128}),
          length(x), pointer(x))
  else
    Base.norm(x)
  end
end

@doc """ sum of vector x """
function sum(
    x :: Vector
)
  if use_SPMP
    ccall((:sum, LIB_PATH), Cdouble,
          (Cint, Ptr{Cdouble}),
          length(x), x)
  else
    Base.sum(x)
  end
end

function minimum(
    x :: Vector
)
  if use_SPMP
    ccall((:minimum, LIB_PATH), Cdouble,
          (Cint, Ptr{Cdouble}),
          length(x), x)
  else
    Base.minimum(x)
  end
end

function abs!(
    w :: Vector,
    x :: Vector
)
  if use_SPMP
    ccall((:CSR_abs, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), w, x)
    w
  else
    w[:] = Base.abs(x)
  end
end

function abs!(
    w :: Vector{Float64},
    x :: Vector{Complex128}
)
  if use_SPMP
    ccall((:CSR_abs_complex, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Complex128}),
          length(x), w, x)
  else
    w[:] = Base.abs(x)
  end
end

function abs!(
    w :: Vector{Float64},
    x :: Vector{Complex128}
)
  if use_SPMP
    ccall((:CSR_abs_complex, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Complex128}),
          length(x), w, x)
  else
    w[:] = abs(x)
  end
end

function exp!(
    w :: Vector,
    x :: Vector
)
  if use_SPMP
    ccall((:CSR_exp, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), w, x)
    w
  else
    w[:] = Base.exp(x)
  end
end

function log1p!(
    w :: Vector,
    x :: Vector
)
  if use_SPMP
    ccall((:CSR_log1p, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), w, x)
    w
  else
    w[:] = Base.log1p(x)
  end
end

function min!(
    w     :: Vector,
    x     :: Vector,
    alpha :: Number
)
  assert(length(x) == length(w))

  if use_SPMP
    ccall((:min, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble),
          length(x), w, x, alpha)
    w
  else
    w[:] = Base.min(x, alpha)
  end
end

@doc """ copy x to y """
function copy!(
    y :: Vector,
    x :: Vector
)
  WAXPB!(y, 1, x, 0)
end

@doc """ w = x./y """
function element_wise_divide!(
    w :: Vector, 
    x :: Vector, 
    y :: Vector
)
  assert(length(x) == length(y))

  if use_SPMP
    ccall((:pointwiseDivide, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(w), pointer(x), pointer(y))
    w
  else
    w[:] = x./y
  end
end

@doc """" Alocate space for and return element-wise division of two vectors. """
function element_wise_divide(
    x :: Vector, 
    y :: Vector
)
  w = Array(Cdouble, length(x))
  element_wise_divide!(w, x, y)
  w
end

@doc """" w = x.*y """
function element_wise_multiply!(
    w :: Vector, 
    x :: Vector, 
    y :: Vector
)
  assert(length(x) == length(y))

  if use_SPMP
    ccall((:pointwiseMultiply, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(w), pointer(x), pointer(y))
    w
  else
    w[:] = x.*y
  end
end

@doc """" w = x.*y """
function element_wise_multiply!(
    w :: Vector{Complex128}, 
    x :: Vector{Float64}, 
    y :: Vector{Complex128}
)
  assert(length(x) == length(y))

  if use_SPMP
    ccall((:pointwiseMultiplyRealWithComplex, LIB_PATH), Void,
          (Cint, Ptr{Complex128}, Ptr{Float64}, Ptr{Complex128}),
          length(x), pointer(w), pointer(x), pointer(y))
  else
    w[:] = x.*y
  end
end

@doc """" w = x.*y """
function element_wise_multiply!(
    w :: Vector{Complex128}, 
    x :: Vector{Float64}, 
    y :: Vector{Complex128}
)
  assert(length(x) == length(y))

  if use_SPMP
    ccall((:pointwiseMultiplyRealWithComplex, LIB_PATH), Void,
          (Cint, Ptr{Complex128}, Ptr{Float64}, Ptr{Complex128}),
          length(x), pointer(w), pointer(x), pointer(y))
  else
    w[:] = x.*y
  end
end

@doc """" Alocate space for and return element-wise multiplication of two vectors. """
function element_wise_multiply(
    x :: Vector{Float64}, 
    y :: Vector{Float64}
)
  w = Array(Cdouble, length(x))
  element_wise_multiply!(w, x, y)
  w
end

@doc """" Alocate space for and return element-wise multiplication of two vectors. """
function element_wise_multiply(
    x :: Vector{Float64}, 
    y :: Vector{Complex128}
)
  w = Array(Complex128, length(x))
  element_wise_multiply!(w, x, y)
  w
end


@doc """ w = alpha*x + beta*y """
function WAXPBY!(
    w     :: Vector, 
    alpha :: Number, 
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector
)
  assert(length(x) == length(y))
  assert(length(x) == length(w))

  if use_SPMP
    ccall((:waxpby, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
          length(x), pointer(w), alpha, pointer(x), beta, pointer(y))
    w
  else
    w[:] = alpha*x + beta*y
  end
end

@doc """ w = alpha*x + beta*y """
function WAXPBY!(
    w     :: Vector{Complex128}, 
    alpha :: Number, 
    x     :: Vector{Complex128}, 
    beta  :: Number, 
    y     :: Vector{Complex128}
)
  assert(length(x) == length(y))
  assert(length(x) == length(w))

  if use_SPMP
    ccall((:waxpby_complex, LIB_PATH), Void,
          (Cint, Ptr{Complex128}, Complex128, Ptr{Complex128}, Complex128, Ptr{Complex128}),
          length(x), pointer(w), alpha, pointer(x), beta, pointer(y))
  else
    w[:] = alpha*x + beta*y
  end
end

@doc """ w = alpha*x + beta*y """
function WAXPBY!(
    w     :: Vector{Complex128}, 
    alpha :: Number, 
    x     :: Vector{Complex128}, 
    beta  :: Number, 
    y     :: Vector{Complex128}
)
  assert(length(x) == length(y))
  assert(length(x) == length(w))

  if use_SPMP
    ccall((:waxpby_complex, LIB_PATH), Void,
          (Cint, Ptr{Complex128}, Complex128, Ptr{Complex128}, Complex128, Ptr{Complex128}),
          length(x), pointer(w), alpha, pointer(x), beta, pointer(y))
  else
    w[:] = alpha*x + beta*y
  end
end

@doc """ Allocate space for and return alpha*x + beta*y """
function WAXPBY(
    alpha :: Number,
    x     :: Vector{Float64}, 
    beta  :: Number, 
    y     :: Vector{Float64}
)
  w = Array{Cdouble}(length(x))
  WAXPBY!(w, alpha, x, beta, y)
  w
end

@doc """ Allocate space for and return alpha*x + beta*y """
function WAXPBY(
    alpha :: Number,
    x     :: Vector{Complex128}, 
    beta  :: Number, 
    y     :: Vector{Complex128}
)
  w = Array(Complex128, length(x))
  WAXPBY!(w, alpha, x, beta, y)
  w
end

@doc """ Allocate space for and return alpha*x + beta*y """
function WAXPBY(
    alpha :: Number,
    x     :: Vector{Complex128}, 
    beta  :: Number, 
    y     :: Vector{Complex128}
)
  w = Array(Complex128, length(x))
  WAXPBY!(w, alpha, x, beta, y)
  w
end

@doc """ w = alpha*x + beta """
function WAXPB!(
    w     :: Vector, 
    alpha :: Number, 
    x     :: Vector, 
    beta  :: Number
)
  assert(length(w) == length(x))

  if use_SPMP
    ccall((:waxpb, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
          length(x), pointer(w), alpha, pointer(x), beta)
    w
  else
    w[:] = alpha*x + beta
  end
end

@doc """ w = alpha*x + beta """
function WAXPB!(
    w     :: Vector{Complex128}, 
    alpha :: Number, 
    x     :: Vector{Complex128}, 
    beta  :: Number
)
  assert(length(w) == length(x))

  if use_SPMP
    ccall((:waxpb_complex, LIB_PATH), Void,
          (Cint, Ptr{Complex128}, Complex128, Ptr{Complex128}, Complex128),
          length(x), pointer(w), alpha, pointer(x), beta)
  else
    w[:] = alpha*x + beta
  end
end

@doc """ w = alpha*x + beta """
function WAXPB!(
    w     :: Vector{Complex128}, 
    alpha :: Number, 
    x     :: Vector{Complex128}, 
    beta  :: Number
)
  assert(length(w) == length(x))

  if use_SPMP
    ccall((:waxpb_complex, LIB_PATH), Void,
          (Cint, Ptr{Complex128}, Complex128, Ptr{Complex128}, Complex128),
          length(x), pointer(w), alpha, pointer(x), beta)
  else
    w[:] = alpha*x + beta
  end
end

@doc """ Allocate space for and return alpha*x + beta """
function WAXPB(
    alpha :: Number, 
    x     :: Vector, 
    beta  :: Number
)
  w = Array(Cdouble, length(x))
  WAXPB!(w, alpha, x, beta)
  w
end
 
@doc """
Compute A * D * B, where A and B have constant structures, D is a diagonal
matrix.

The A and B are inputed as AT and BT (transposed A and B) in Julia CSC format,
which become A and B in CSR.

An additional fknob is passed in, which contains a matrix knob representing the
output of this function, and 3 other matrix knobs representing the input of this
function (A, D, B).

In addition to its original semantics, which is returning A * D *B, this
function also reuses context info (from the input matrix knobs), and generates
new context info (for the output matrix knob). 

If the source statements are 
    for  
        X = A * D * B
they have been transformed into the following form before executing this function:
    C = Symbol("ouptut") # representing the ouptut of the call at this call site
    mknob_C = new_matrix_knob(properties of C)
    mknob_A = new_matrix_knob(properties of A)
    mknob_D = new_matrix_knob(properties of D)
    mknob_B = new_matrix_knob(properties of B)
    mknob_X = new_matrix_knob(properties of X)
    fknob_ADB = new_function_knob()
    add mknob_C, _A, _D, _B to fknob_ADB
    for 
        X = ADB(AT, D, BT, fknob_ADB)
        knob_X = mknob_C  #propagate mknob_C information to mknob_X.
                          #This may or may not be one to one copy. For example,
                          # if C is constant valued, but X is only constant
                          # structured, then we may not pass C's memory
                          # to X directly; instead, we might pass a copy of C's
                          # memory.

In general, in context-sensitive optimizations, there are two problems:
(1) how a function generates and reuses context info?
(2) how to propagate context info? 
Here we focus on the first problem.

ADB(A', D, B', fknob_ADB) functionality:
    if fknob_ADB == NULL
        return A * D * B
    
    get the output mknob_C, and input mknob_A, D, B from fknob_ADB
    assert(mknob_C != NULL)
    assert(mknob_A != NULL)
    assert(mknob_D != NULL)
    assert(mknob_B != NULL)
    
    # Check inputs
    assert(mknob_A->constant_structured)
    assert(mknob_D->diagonal)
    assert(mknob_B->constant_structured)
    
    # Compute output
    if mknob_C->matrix == NULL
        mknob_C->matrix = adb_inspect(A', B')
    mknob_C->matrix = CSR_ADB(mknob_C->matrix,  A', B', diag(D))

    # Return the result.
    return mknob_C->matrix in CSC format
"""
function ADB(
    A     :: SparseMatrixCSC,
    D     :: SparseMatrixCSC,
    B     :: SparseMatrixCSC,
    fknob :: Ptr{Void}
)
    if fknob == C_NULL
        # This case should never be reached. We leave it here only to show how
        # a library function can have the same interface with and without 
        # context info.
        assert(false)
        return AT' * D * BT'
    end

    mknob_C = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    mknob_A = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 1)
    mknob_D = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 2)
    mknob_B = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 3)
    assert(mknob_C != C_NULL)
    assert(mknob_A != C_NULL)
    assert(mknob_D != C_NULL)
    assert(mknob_B != C_NULL)

    # Check inputs and output
    assert(ccall((:IsConstantStructured, LIB_PATH), Bool, (Ptr{Void},), mknob_A))
    #TODO: assert(mknob_D->diagonal)
    assert(ccall((:IsConstantStructured, LIB_PATH), Bool, (Ptr{Void},), mknob_B))

    fill_matrix_info_once(mknob_A, A)
    fill_matrix_info_once(mknob_D, D)
    fill_matrix_info_once(mknob_B, B)

    m = size(A, 1)
    n = size(B, 2)
    k = size(A, 2)

    C_colptr_ref = Ref{Ptr{Cint}}(C_NULL)
    C_rowval_ref = Ref{Ptr{Cint}}(C_NULL)
    C_nzval_ref = Ref{Ptr{Cdouble}}(C_NULL)

    ccall((:ADB, LIB_PATH), Void,
         (Cint, Cint, Cint,
           Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, Ref{Ptr{Cdouble}},
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
           Ptr{Cdouble},
           Ptr{Void}),
         m, n, k,
         C_colptr_ref, C_rowval_ref, C_nzval_ref,
         A.colptr, A.rowval, A.nzval,
         B.colptr, B.rowval, B.nzval,
         D.nzval,
         fknob)

    C_colptr = pointer_to_array(C_colptr_ref[], n + 1)
    nnz = C_colptr[n + 1] - 1
    C_rowval = pointer_to_array(C_rowval_ref[], nnz)
    C_nzval = pointer_to_array(C_nzval_ref[], nnz)
    SparseMatrixCSC{Cdouble, Cint}(m, n, C_colptr, C_rowval, C_nzval)
end

# A*A'
function SpSquareWithEps(
    A     :: SparseMatrixCSC{Float64, Int32},
    eps   :: Number,
    fknob :: Ptr{Void} = C_NULL)

    if fknob != C_NULL
        mknob_A = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
        assert(mknob_A != C_NULL)
        fill_matrix_info_once(mknob_A, A)
    end
    
    m = size(A, 1)
    n = size(A, 2)

    C_colptr_ref = Ref{Ptr{Cint}}(C_NULL)
    C_rowval_ref = Ref{Ptr{Cint}}(C_NULL)
    C_nzval_ref = Ref{Ptr{Cdouble}}(C_NULL)

    ccall((:SpSquareWithEps, LIB_PATH), Ptr{Void},
          (Cint, Cint,
           Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, Ref{Ptr{Cdouble}},
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
           Cdouble,
           Ptr{Void}),
          m, n,
          C_colptr_ref, C_rowval_ref, C_nzval_ref,
          A.colptr, A.rowval, A.nzval,
          eps,
          fknob)

    C_colptr = pointer_to_array(C_colptr_ref[], n + 1)
    nnz = C_colptr[n + 1] - 1
    C_rowval = pointer_to_array(C_rowval_ref[], nnz)
    C_nzval = pointer_to_array(C_nzval_ref[], nnz)
    SparseMatrixCSC{Cdouble, Cint}(m, n, C_colptr, C_rowval, C_nzval)
end

# alpha*A + beta*B
function SpAdd(
    alpha :: Number,
    A     :: SparseMatrixCSC{Float64, Int32},
    beta  :: Number,
    B     :: SparseMatrixCSC{Float64, Int32})

    m = size(A, 1)
    n = size(A, 2)
    assert(m == size(B, 1))
    assert(n == size(B, 2))

    C_colptr_ref = Ref{Ptr{Cint}}(C_NULL)
    C_rowval_ref = Ref{Ptr{Cint}}(C_NULL)
    C_nzval_ref = Ref{Ptr{Cdouble}}(C_NULL)

    ccall((:SpAdd, LIB_PATH), Ptr{Void},
          (Cint, Cint,
           Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, Ref{Ptr{Cdouble}},
           Cdouble,
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
           Cdouble,
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
           Ptr{Void}),
          m, n,
          C_colptr_ref, C_rowval_ref, C_nzval_ref,
          alpha,
          A.colptr, A.rowval, A.nzval,
          beta,
          B.colptr, B.rowval, B.nzval,
          C_NULL)

    C_colptr = pointer_to_array(C_colptr_ref[], n + 1)
    nnz = C_colptr[n + 1] - 1
    C_rowval = pointer_to_array(C_rowval_ref[], nnz)
    C_nzval = pointer_to_array(C_nzval_ref[], nnz)
    SparseMatrixCSC{Cdouble, Cint}(m, n, C_colptr, C_rowval, C_nzval)
end

function trace(
    A   :: SparseMatrixCSC{Float64, Int32})

    ccall((:Trace, LIB_PATH), Cdouble,
          (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
          size(A, 1), A.colptr, A.rowval, A.nzval)
end

@doc """ 
Compute the sparse Cholesky factorization of a sparse matrix A, where A is
constant in structure.

Functionality:
    if fknob == NULL
        return cholfact_int32(A)
    
    get the output/input mknob_O and A from fknob
    
    # Check input
    assert(mknob_A->constant_structured)
    if mknob_A->matrix == NULL
        return cholfact_int32(A)
        
    if mknob_O->dss_handle == NULL
        mknob_O->dss_handle = dss_analyze(mknob_A->matrix)
    dss_factor(mknob_O->dss_handle, mknob_A->matrix)
"""
function cholfact_int32(
    A     :: SparseMatrixCSC,
    fknob :: Ptr{Void}
)
    if fknob == C_NULL
        # This case should never happen
        assert(false)
        return cholfact_int32(A)
    end

    mknob_A = ccall((:GetMatrixKnob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    assert(mknob_A != C_NULL)

    fill_matrix_info_once(mknob_A, A)

    ccall((:CholFact, LIB_PATH), Ptr{Void},
         (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Void}),
         size(A, 1), size(A, 2), A.colptr, A.rowval, A.nzval, fknob)
end

@doc """ 
Compute y = R \ t, where R is constant in structure.

Functionality:
    if fknob == NULL
       return y = R \ t
   get the input mknob_R from fknob
   if mknob_R->dss_handle == NULL
       return y = R \ t2
   else
       opt = MKL_DSS_DEFAULTS
       dss_solve!(mknob_R->dss_handle, t, y)
"""
function cholfact_inverse_divide!(
    y     :: Vector,
    R     :: Any,
    b     :: Vector,
    fknob :: Ptr{Void}
)
    if fknob == C_NULL
        # This case should never happen
        assert(false)
        y[:] = R \ b
        return y
    end

    ccall((:CholFactInverseDivide, LIB_PATH), Void,
          (Ptr{Void}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
          R, y, b, fknob)
    y
end

@doc """ 
Compute y = R \ t, where R is constant in structure.

Functionality:
    if fknob == NULL
       return y = R \ t
   get the input mknob_R from fknob
   if mknob_R->dss_handle == NULL
       return y = R \ t2
   else
       opt = MKL_DSS_DEFAULTS
       dss_solve!(mknob_R->dss_handle, t, y)
"""
function cholfact_inverse_divide(
    R     :: Any,
    b     :: Vector,
    fknob :: Ptr{Void}
)
    y = zeros(length(b))
    cholfact_inverse_divide!(y, R, b, fknob)
    y
end

# Some symbolic names for each permutation vector.
const NOT_PERM_YET = 0 # Maybe permute, but the color is not decided yet.
const ROW_PERM     = 1
const ROW_INV_PERM = 2
const COL_PERM     = 3
const COL_INV_PERM = 4
const NEVER_PERM   = 5 # Never permute. This state can not change to any other.

const LOG_REORDERING = false

@doc """
Reorder arrays based on the the given permutation and inverse permutation vectors.
If inverse_reorder is true, do inverse reordering instead.
"""
function reorder_arrays(
    inverse_reorder     :: Bool,
    permutation_vectors :: Vector,
    arrays...
)
    if LOG_REORDERING
      println("reorder_arrays")
    end

    matrices_done = false
    array_tuple   = arrays[1]
    i = 1
    while i <= length(array_tuple)
        A = array_tuple[i]
        if A == :__delimitor__
            # A delimiter that separate matrices from vectors
            # It means that all matrices have been processed
            matrices_done = true
            i             = i + 1
            continue
        end
        if matrices_done
            assert(typeof(A) <: AbstractVector)
            # For each vector, 1 permutation vector follows
            perm_vec_idx = array_tuple[i + 1]
            perm         = permutation_vectors[perm_vec_idx]
            i    = i + 2
            if inverse_reorder
                inverse_reorder_vector!(A, perm)
            else
                # reordering(v, ROW_PERM) is equivalent to 
                # reverse_reordering(v, ROW_INV_PERM) but is slower. So change
                # to the second form if possible.
                if perm_vec_idx == ROW_PERM || perm_vec_idx == COL_PERM
                    perm = permutation_vectors[perm_vec_idx + 1]
                    inverse_reorder_vector!(A, perm)
                else
                    assert(perm_vec_idx == ROW_INV_PERM || perm_vec_idx == COL_INV_PERM)
                    perm = permutation_vectors[perm_vec_idx - 1]
                    inverse_reorder_vector!(A, perm)
                end
            end
        else
            assert(typeof(A) <: AbstractSparseMatrix)
            # For each matrix, 2 permutation vectors and 1 matrix knob follow
            perm1 = permutation_vectors[array_tuple[i + 1]]
            perm2 = permutation_vectors[array_tuple[i + 2]]
            mknob = array_tuple[i + 3]
            fill_matrix_info_once(mknob, A)
            i     = i + 4
            if inverse_reorder
                reorder_matrix!(A, perm2, perm1)
            else
                reorder_matrix!(A, perm1, perm2)
            end
        end
    end
end

@doc """
The first time an reordering decision maker (a function, represented by the
fknob) has decided that reordering should be done, this function should be 
called to reorder the specified arrays.

Status is a vector of the following elements:
(1) reordering_done:       Bool, true if reordering has already been done
(2) row_perm:              row permutation vector
(3) row_inv_perm:          row inverse permutaiton vector
(4) col_perm:              column permutation vector
(5) col_inv_perm:          column inverse permutaiton vector
(6) reordering_time:       time spent on reordering so far (for statistics)
The vector is updated when the function returns.
"""
function reordering(
    fknob  :: Ptr{Void},
    status :: Vector,
    arrays...
)
    reordering_done = status[1]
    if !reordering_done
        row_perm, row_inv_perm, col_perm, col_inv_perm = get_reordering_vectors(fknob)
        if length(row_perm) > 0
            assert(length(row_inv_perm) > 0 && length(col_perm) > 0 && length(col_inv_perm) > 0)

            status[2] = row_perm
            status[3] = row_inv_perm
            status[4] = col_perm
            status[5] = col_inv_perm

            reordering_time  = status[6]
            reordering_time -= time()

            reorder_arrays(false, status[2 : 5], arrays)

            reordering_time += time()
            status[6]        = reordering_time
            status[1]        = true # reordering is done 
        end
    end
end

@doc """
Reverse reorder the specified arrays. 

Status is a vector of the following 4 elements. See the above reordering()
function for details.
"""
function reverse_reordering(
    status :: Vector,
    arrays...
)
    reordering_done = status[1]
    if reordering_done
        reordering_time  = status[6]
        reordering_time -= time()

        reorder_arrays(true, status[2 : 5], arrays)

        reordering_time += time()
        status[6]        = reordering_time
    end
end

@doc """
Read a matrix market file. Force it to be symmetric if force_symmetric is true.
Error when symmetric_only is true but the matrix is not symmetric (and not
forced to be symmetric). It calls SpMP library to do the actual work, which is 
faster than a pure Julia version.
"""
function matrix_market_read(
    filename        :: AbstractString,
    symmetric_only  :: Bool = false,
    force_symmetric :: Bool = false
)
    # Read header to figure out if it's an array or sparse matrix
    mmfile = open(filename,"r")
    # Read first line
    firstline = chomp(readline(mmfile))
    tokens = split(firstline)
    if length(tokens) != 5
        throw(ParseError(string("Not enough words on first line: ", ll)))
    end
    if tokens[1] != "%%MatrixMarket"
        throw(ParseError(string("Not a valid MatrixMarket header:", ll)))
    end
    (head1, rep, field, symm) = map(lowercase, tokens[2:5])
    if head1 != "matrix"
        throw(ParseError("Unknown MatrixMarket data type: $head1 (only \"matrix\" is supported)"))
    end
    close(mmfile)

    if rep == "array"
      m_ref = Ref{Cint}(0)
      n_ref = Ref{Cint}(0)
      v_ref = Ref{Ptr{Cdouble}}(C_NULL)

      ccall((:load_vector_matrix_market, LIB_PATH), Void,
           (Ptr{Cchar}, Ref{Ptr{Cdouble}}, Ref{Cint}, Ref{Cint}),
           filename, v_ref, m_ref, n_ref)
      assert(n_ref[] == 1)

      return pointer_to_array(v_ref[], m_ref[])
    else
      # Read into a COO array
      sizes::Vector{Cint} = [0, 0, 0, 0]
      ccall((:load_matrix_market_step1, LIB_PATH), Void, 
          (Ptr{UInt8}, Ptr{Cint}, Bool, Bool), filename, pointer(sizes), force_symmetric, false)

      is_symmetric::Cint = sizes[1] # 0/1: true/false    
      if symmetric_only && is_symmetric == 0
          # We cannot handle asymmetric matrix right now, since we need a matrix
          # to be symmetric, and thus we can use its CSC representation as CSR (The
          # SpMP library is based on CSR, while Julia is based on CSC.)
          throw(AsymmetricMatrixMarketFile(filename))
      end

      n::Cint = sizes[2]
      m::Cint = sizes[3]
      assert(is_symmetric == 0 || m == n)
      nnz::Cint = sizes[4]  # Note: if symmetric, nnz includes elements in both 
                            # upper and lower triangle 

      v = Array(Cdouble, nnz)
      i = Array(Cint, nnz)
      j = Array(Cint, n + 1)
      
      A = SparseMatrixCSC{Cdouble, Cint}(m, n, j, i, v)

      # Convert the COO array to CSC array.
      ccall((:load_matrix_market_step2, LIB_PATH), Void, 
          (Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Bool),
          filename, pointer(j), pointer(i), pointer(v), pointer(sizes), true)

      return A
    end
end

function lbfgs_compute_direction!(
    dx,
    k     :: Int,
    it    :: Int,
    n     :: Int,
    S,
    Y,
    dfk)

    if use_SPMP
      ccall((:LBFGSComputeDirection, LIB_PATH), Void,
            (Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            k, it, n, S, Y, dfk, dx)
    else
      a = zeros(k)

      mm = min(it - 1, k)
      begin_idx = it + k - mm - 2

      yk = 1
      if it <= k
      else
        skm1 = S[:,(begin_idx + mm)%k + 1]
        ykm1 = Y[:,(begin_idx + mm)%k + 1]
        yk = dot(skm1, ykm1)/dot(ykm1, ykm1)
      end

      p = zeros(mm, 1)
      for i = 1:mm
        p[i] = 1.0/dot(S[:,(begin_idx + i)%k + 1], Y[:,(begin_idx + i)%k + 1])
      end

      q = copy(dfk)
      for i = mm:-1:1
        si = S[:,(begin_idx + i)%k + 1]
        yi = Y[:,(begin_idx + i)%k + 1]
        a[i] = p[i]*dot(si, q)
        q = q - a[i]*yi
      end

      r = yk*q
      for i = 1:mm
        si = S[:,(begin_idx + i)%k + 1]
        yi = Y[:,(begin_idx + i)%k + 1]
        b = p[i]*dot(yi, r)
        r = r + si*(a[i] - b)
      end

      dx[:] = -r
    end
end 

function lbfgs_loss_function1(yXw::Vector, w::Vector, lambda::Number)
  m = length(yXw)

  if use_SPMP
    ccall((:LBFGSLossFunction1, LIB_PATH), Cdouble,
          (Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble),
          m, length(w), yXw, w, lambda)
  else
    s = sum(log1p(exp(-abs(yXw))) - min(yXw, 0))
    s/m + (lambda/2)*norm(w)^2
  end 
end

function lbfgs_loss_function2!(temp::Vector, y::Vector, yXw::Vector)
  if use_SPMP
    ccall((:LBFGSLossFunction2, LIB_PATH), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          temp, length(y), y, yXw)
  else
    temp[:] = y./(1 + exp(yXw))
  end
end

function get_spmp_spmv_time()
  ccall((:GetSpMPSpMVTime, LIB_PATH), Cdouble, ())
end

function reset_spmp_spmv_time()
  ccall((:ResetSpMPSpMVTime, LIB_PATH), Void, ())
end

function get_spmp_trisolve_time()
  ccall((:GetSpMPTriangularSolveTime, LIB_PATH), Cdouble, ())
end

function reset_spmp_trisolve_time()
  ccall((:ResetSpMPTriangularSolveTime, LIB_PATH), Void, ())
end

function get_knob_spmv_time()
  ccall((:GetKnobSpMVTime, LIB_PATH), Cdouble, ())
end

function reset_knob_spmv_time()
  ccall((:ResetKnobSpMVTime, LIB_PATH), Void, ())
end

function get_knob_trisolve_time()
  ccall((:GetKnobTriangularSolveTime, LIB_PATH), Cdouble, ())
end

function reset_knob_trisolve_time()
  ccall((:ResetKnobTriangularSolveTime, LIB_PATH), Void, ())
end

function set_knob_log_level(level)
  ccall((:SetLogLevel, LIB_PATH), Void, (Cint,), level)
end

function set_reuse_inspection(reuse)
  ccall((:SetReuseInspection, LIB_PATH), Void, (Bool,), reuse)
end

function set_collective_structure_prediction(collective)
  ccall((:SetCollectiveStructurePrediction, LIB_PATH), Void, (Bool,), collective)
end
