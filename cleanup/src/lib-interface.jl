# This file contains the routines that interfaces with SpMP library.

const LIB_PATH = libcsr

@doc """ Create a new matrix knob. """
function new_matrix_knob()
    mknob = ccall((:NewMatrixKnob, LIB_PATH), Ptr{Void}, ())
end

@doc """ Increment the version of a matrix knob. """
function increment_matrix_version(
    mknob :: Ptr{Void}
)
    ccall((:IncrementMatrixVersion, LIB_PATH), Void, (Ptr{Void},), mknob)
end

@doc """ Delete a matrix knob. """
function delete_matrix_knob(
    mknob :: Ptr{Void}
)
    ccall((:DeleteMatrixKnob, LIB_PATH), Void, (Ptr{Void},), mknob)
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
    fknob_creator :: String
)
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

@doc """ Delete a function knob. """
function delete_function_knob(
    fknob_deletor :: String, 
    fknob         :: Ptr{Void}
)
    # Simulate ccall((fknob_deletor, LIB_PATH), Ptr{Void}, ()). We cannot directly
    # use this ccall here because (fknob_deletor, LIB_PATH) is treated as a tuple,
    # instead of a pointer or expression.
    expr = Expr(:call, 
        TopNode(:ccall), 
        Expr(:call, TopNode(:tuple), QuoteNode(fknob_deletor), LIB_PATH), 
        Void, 
        Expr(:call, TopNode(:svec))
    )
    eval(expr)
end

@doc """ Context-sensitive forward triangular solver. """
function fwdTriSolve!(
    A     :: SparseMatrixCSC, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    y = zeros(Cdouble, length(b))
    ccall((:ForwardTriangularSolve, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(y), pointer(b), fknob)
    b = copy(y)
end

@doc """ Context-sensitive backward triangular solver. """
function bwdTriSolve!(
    A     :: SparseMatrixCSC, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    y = zeros(Cdouble, length(b))
    ccall((:BackwardTriangularSolve, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(y), pointer(b), fknob)
    b = copy(y)
end

@doc """ 
Reorder a sparse matrix A and store the result in new_A. A itself is not changed.
If get_permutation is true, then compute the permutation and inverse permutation
vector P and inverse_P. Otherwise, P and inverse_P are given inputs.
One_based_input/output tells if A and new_A are 1-based.
"""
function reorder_matrix(
    A                :: SparseMatrixCSC, 
    new_A            :: SparseMatrixCSC, 
    P                :: Vector, 
    inverse_P        :: Vector, 
    get_permutation  :: Bool, 
    one_based_input  :: Bool, 
    one_based_output :: Bool
)
    ccall((:CSR_ReorderMatrix, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool, Bool, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(new_A.colptr), pointer(new_A.rowval), pointer(new_A.nzval),
               pointer(P), pointer(inverse_P), get_permutation, one_based_input, 
               one_based_output)
end

@doc """ Reorder vector V into new vector new_V with permutation vector P """
function reorder_vector(
    V     :: Vector, 
    new_V :: Vector, 
    P     :: Vector
)
    ccall((:reorderVector, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(new_V), pointer(P), length(V))
end

@doc """ Reversely reorder vector V into new vector new_V with permutation vector P """
function reverse_reorder_vector(
    V     :: Vector, 
    new_V :: Vector, 
    P     :: Vector
)
    ccall((:reorderVectorWithInversePerm, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(new_V), pointer(P), length(V))
end

@doc """
Return a new representation of the sparse matrix, which is in CSC 
format, in CSR format.
"""
function create_CSR(A :: SparseMatrixCSC)
    ccall((:CSR_Create, LIB_PATH), Ptr{Void},
       (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint),
       A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval), 1)
end

@doc """ Destroy the CSR representation of the sparse matrix. """
function destroy_CSR(A :: Ptr{Void})
    ccall((:CSR_Destroy, LIB_PATH), Void, (Ptr{Void},), A)
end

@doc """ w = alpha*A*x + beta*y + gamma """
function SpMV!(
    w     :: Vector, 
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector, 
    gamma :: Number
)
    assert(length(w) == length(y))

    if use_SPMP
        A1 = create_CSR(A)
        ccall((:CSR_MultiplyWithVector, LIB_PATH), Void,
              (Ptr{Cdouble}, Cdouble, Ptr{Void}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
              pointer(w), alpha, A1, pointer(x), beta, pointer(y), gamma)
        destroy_CSR(A1)
    else
        # use Julia implementation
        w = alpha * A * x + beta * y + gamma
    end
    w
end

@doc """ y = A*x """
SpMV!(y :: Vector, A :: SparseMatrixCSC, x :: Vector) = 
    SpMV!(y, one(eltype(x)), A, x, zero(eltype(y)), y, zero(eltype(x)))

@doc """ alpha*A*x + beta*y + gamma """
function SpMV(
    alpha :: Number, 
    A     :: SparseMatrixCSC, 
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector, 
    gamma :: Number
)
   w = Array(Cdouble, length(x))
   SpMV!(w, alpha, A, x, beta, y, gamma)
   w
end

@doc """ alpha*A*x + beta*y """
SpMV(alpha :: Number, A :: SparseMatrixCSC, x :: Vector, beta :: Number, y :: Vector) =
    SpMV(alpha, A, x, beta, y, zero(eltype(x)))

@doc """ alpha*A*x + y """
SpMV(alpha :: Number, A :: SparseMatrixCSC, x :: Vector, y :: Vector) = 
    SpMV(alpha, A, x, one(eltype(y)), y, zero(eltype(x)))

@doc """ alpha*A*x """
SpMV(alpha :: Number, A :: SparseMatrixCSC, x :: Vector) = 
    SpMV(alpha, A, x, zero(eltype(x)), x, zero(eltype(x)))

@doc """ A*x """
SpMV(A :: SparseMatrixCSC, x :: Vector) = 
    SpMV(one(eltype(x)), A, x, zero(eltype(x)), x, zero(eltype(x)))

@doc """
SpMV. In addition, measure the potential benefit of reordering its input matrix,
if reordering has not done yet.
"""
function SpMV_conditional_reordering(
    A                  :: SparseMatrixCSC, 
    x                  :: Vector, 
    first_reorder_done :: Bool, 
    beneficial         :: Vector{Bool}
)
    if first_reorder_done
        return SpMV(A, x)
    end
    time1 = time()
    y = SpMV(A, x)
    time2 = time()
    nnz = size(A.nzval, 1)
    rows = A.m
    SpMV_bandwidth = (nnz * 12 + rows * 3 * 8) / (time2 - time1) * 128 / 10 ^ 9
    # ISSUE: here the machine peak bandwidth and the distance are fixed constants.
    # The machine peak bandwidth should be set after running STREAM benchmark once
    # during the installation time of the SparseAccelerator. The distance should
    # be changed to make it portable for different machines. 
    machine_peak_bandwidth = 60
    if abs(SpMV_bandwidth - machine_peak_bandwidth) > 15
        beneficial[1] = true
    end
    return y
end

@doc """ Initialize the conditions for conditional reordering. """
function init_conditional_reordering(
    first_reorder_done :: Bool, 
    beneficial         :: Vector
)
    # Benefical can be just a scalar var. But we will treat it as a 1-element
    # array so that we do not have difficulty in changing it in calling 
    # SpMV_conditional_reordering
    beneficial = Vector{Bool}()
    push!(beneficial, false)
    first_reorder_done = false
end

@doc """ Dot product of vector x and y """
function dot(
    x :: Vector, 
    y :: Vector
)
  assert(length(x) == length(y))

  if use_SPMP
    ccall((:dot, LIB_PATH), Cdouble,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(x), pointer(y))
  else
    dot(x, y)
  end
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
  else
    w = x./y
  end
  w
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
  else
    w = x.*y
  end
  w
end

@doc """" Alocate space for and return element-wise multiplication of two vectors. """
function element_wise_multiply(
    x :: Vector, 
    y :: Vector
)
  w = Array(Cdouble, length(x))
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
  else
    w = alpha*x + beta*y
  end
  w
end

@doc """ Allocate space for and return alpha*x + beta*y """
function WAXPBY(
    alpha :: Number,
    x     :: Vector, 
    beta  :: Number, 
    y     :: Vector
)
  w = Array(Cdouble, length(x))
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
  else
    w = alpha*x + beta
  end
  w
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
Read a matrix market file. Force it to be symmetric if force_symmetric is true.
It calls SpMP library to do the actual work, which is faster than a pure
Julia version.
"""
function matrix_market_read(
    filename        :: String, 
    force_symmetric :: Bool = false
)
    # Read into a COO array
    sizes::Vector{Cint} = [0, 0, 0, 0]
    ccall((:load_matrix_market_step1, LIB_PATH), Void, 
        (Ptr{Uint8}, Ptr{Cint}, Bool), filename, pointer(sizes), force_symmetric)

    is_symmetric::Cint = sizes[1] # 0/1: true/false    
    if is_symmetric == 0
        # We cannot handle asymmetric matrix right now, since we need a matrix
        # to be symmetric, and thus we can use its CSC representation as CSR (The
        # SpMP library is based on CSR, while Julia is based on CSC.)
        throw(AsymmetricMatrixMarketFile(filename))
    end

    m::Cint = sizes[2]
    n::Cint = sizes[3]
    assert(m ==n)
    nnz::Cint = sizes[4]  # Note: if symmetric, nnz includes elements in both 
                          # upper and lower triangle 

    v = Array(Cdouble, nnz)
    i = Array(Cint, nnz)
    j = Array(Cint, n + 1)
    
    A = SparseMatrixCSC{Cdouble, Cint}(m, n, j, i, v)

    # Convert the COO array to CSC array.
    ccall((:load_matrix_market_step2, LIB_PATH), Void, 
        (Ptr{Uint8}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Bool),
        filename, pointer(j), pointer(i), pointer(v), pointer(sizes), true)

    A
end