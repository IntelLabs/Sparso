# This file contains the routines that interfaces with SpMP library.

const LIB_PATH = libcsr

@doc """ Create a new matrix knob. """
function new_matrix_knob(
    A     :: SparseMatrixCSC
 )
    mknob = ccall((:NewMatrixKnob, LIB_PATH), Ptr{Void},
                   (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
                   A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval))
end

@doc """ Increment the version of a matrix. """
function increment_matrix_version(
    mknob :: Ptr{Void}
)
    ccall((:IncrementMatrixVersion, LIB_PATH), Void, (Ptr{Void},), mknob)
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

const DERIVATIVE_TYPE_TRANSPOSE = 0
const DERIVATIVE_TYPE_SYMMETRIC = 1
const DERIVATIVE_TYPE_LOWER_TRIANGULAR = 2
const DERIVATIVE_TYPE_UPPER_TRIANGULAR = 3

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
    L     :: SparseMatrixCSC,
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    y = zeros(Cdouble, length(b))
    ccall((:ForwardTriangularSolve, LIB_PATH), Void,
           (
            Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void},
           ),
            L.m, L.n, pointer(L.colptr), pointer(L.rowval), pointer(L.nzval),
            pointer(y), pointer(b), fknob)
#println("\t\tFwdTriSolve! done: sum y = ", sum(y), " y=", y)
#println("\t\tFwdTriSolve! done: sum b = ", sum(b), " b=", b)
    b[:] = y
end

@doc """ 
Context-sensitive backward triangular solver. 
"""
function bwdTriSolve!(
    U     :: SparseMatrixCSC, 
    b     :: Vector,
    fknob :: Ptr{Void}
 )
    y = zeros(Cdouble, length(b))
    ccall((:BackwardTriangularSolve, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
               U.m, U.n, pointer(U.colptr), pointer(U.rowval), pointer(U.nzval),
               pointer(y), pointer(b), fknob) #fknob)
#println("\t\tBwdTriSolve! done: sum y = ", sum(y))#, " y=", y)
#println("\t\tBwdTriSolve! done: sum b = ", sum(b), " b=", b)
    b[:] = y
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
A * D * B, where D is a diagonal matrix, and A, D, and B all have constant structures.
A and B are given in transposed form, as they are in CSC in Julia but libcsr
takes only CSR.

ADB(A', D, B', fknob_ADB) functionality:
   if fknob_ADB == NULL
       return A * D * B
   get the input mknob_A and mknob_B from fknob_ADB
   if !mknob_A->constant_structured || !mknob_B->constant_structured (This check is actually already done by pattern matching)
       return A * D * B
   get the output mknob from fknob_ADB
   if output mknob->structure_proxy == NULL
       ADAT = adb_inspect(A', A), 
       output mknob->structure_proxy = ADAT
       output mknob->matrix = ADAT
   d = diag(D)
   adb_execute!(output mknob->matrix, A', A, d)
"""
function ADB(
    AT    :: SparseMatrixCSC,
    D     :: SparseMatrixCSC,
    BT    :: SparseMatrixCSC,
    fknob :: Ptr{Void}
)
    if fknob == NULL
        return A * D * B
    end
    
    # Skip the checking of the input mknobs as it has been done by pattern matching
    
    mknob_out = ccall((:GetMknob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    structure_proxy = ccall((:GetStructureProxy, LIB_PATH), Ptr{Void}, (Ptr{Void}, ), mknob_out)
    csrA = SparseAccelerator.create_CSR(AT)
    csrB = SparseAccelerator.create_CSR(BT)
    if structure_proxy == C_NULL
        csrADB = ccall((:ADBInspect, LIB_PATH), Ptr{Void},
                        (Ptr{Void}, Ptr{Void}, Ptr{Void}),
                         csrA, csrB)
        ccall((:SetStructureProxy, LIB_PATH), Void, (Ptr{Void}, Ptr{Void}), mknob_out, csrADB)
        ccall((:SetMatrix, LIB_PATH), Void, (Ptr{Void}, Ptr{Void}), mknob_out, csrADB)
    else
        csrADB = ccall((:GetMatrix, LIB_PATH), Ptr{Void}, (Ptr{Void}, ), mknob_out)
    end
    
    d = diag(D)
    ccall((:CSR_ADB, LIB_PATH), Void,
          (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Cdouble}),
           csrADB, csrA, csrB, d)

    # Represent the result in CSC format
    m = size(A, 1)
    rowptr = pointer_to_array(
        ccall((:CSR_GetRowPtr, LIB_PATH), Ptr{Cint}, (Ptr{Void},), csrADB), (m + 1,))
    nnz = rowptr[m + 1] - 1
    colidx = pointer_to_array(
        ccall((:CSR_GetColIdx, LIB_PATH), Ptr{Cint}, (Ptr{Void},), csrADB), (nnz,))
    values = pointer_to_array(
        ccall((:CSR_GetValues, LIB_PATH), Ptr{Cdouble}, (Ptr{Void},), csrADAT), (nnz,))
    ADB = SparseMatrixCSC{Cdouble, Cint}(m, m, rowptr, colidx, values)

    destroy_CSR(csrA)
    destroy_CSR(csrB)

    ADB
end

const MKL_DSS_DEFAULTS = 0
const MKL_DSS_NON_SYMMETRIC = 536871104
const MKL_DSS_SUCCESS = 0
const MKL_DSS_AUTO_ORDER = 268435520
const MKL_DSS_POSITIVE_DEFINITE = 134217792

function dss_analyze(A :: SparseMatrixCSC)
    handle = Int[0]
    opt = MKL_DSS_DEFAULTS
    error = ccall((:dss_create, LIB_PATH), Cint,
        (Ptr{Void}, Ptr{Cint}),
        handle, &opt)
    if error != MKL_DSS_SUCCESS
        println("dss_create returned error code $error")
    end

    opt = MKL_DSS_NON_SYMMETRIC
    m = size(A, 1)
    nnz = A.colptr[m + 1] - 1
    error = ccall((:dss_define_structure, LIB_PATH), Cint,
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        handle, &opt, A.colptr, &m, &m, A.rowval, &nnz)
    if error != MKL_DSS_SUCCESS
        println("dss_define_structure returned error code $error")
    end

    opt = MKL_DSS_AUTO_ORDER
    error = ccall((:dss_reorder, LIB_PATH), Cint,
                (Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
                handle, &opt, C_NULL)
    if error != MKL_DSS_SUCCESS
        println("dss_reorder returned error code $error")
    end

    handle
end

function dss_factor(handle, A::SparseMatrixCSC)
    opt = MKL_DSS_POSITIVE_DEFINITE
    error = ccall((:dss_factor_real, LIB_PATH), Cint,
                (Ptr{Void}, Ptr{Cint}, Ptr{Cdouble}),
                handle, &opt, A.nzval)
    if error != MKL_DSS_SUCCESS
        println("dss_factor_real returned error code $error")
    end
end

# solve A*sol = rhs
function dss_solve!(handle, rhs::Vector, sol::Vector)
    opt = MKL_DSS_DEFAULTS
    nrhs = 1
    error = ccall((:dss_solve_real, LIB_PATH), Cint,
                (Ptr{Void}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
                handle, &opt, rhs, &nrhs, sol)
end

@doc """ 
Compute the sparse Cholesky factorization of a sparse matrix A. 

Functionality:
    if fknob == NULL
        return cholfact_int32(A)
    get the input mknob from fknob
    if input mknob->structure_proxy == NULL || input mknob->matrix == NULL
        return cholfact_int32(B)
    get the output mknob from fknob
    if output mknob->dss_handle == NULL
        dss_handle = dss_analyze(input mknob->structure_proxy)
        output mknob->dss_handle = dss_handle
    dss_factor(output mknob->dss_handle, input mknob->matrix)
"""
function cholfact_int32(
    A     :: SparseMatrixCSC,
    fknob :: Ptr{Void}
)
    if fknob == NULL
        return cholfact_int32(A)
    end

    mknob_in        = ccall((:GetMknob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 1)
    structure_proxy = ccall((:GetStructureProxy, LIB_PATH), Ptr{Void}, (Ptr{Void}, ), mknob_in)
    matrix          = ccall((:GetMatrix, LIB_PATH), Ptr{Void}, (Ptr{Void}, ), mknob_in)
    if structure_proxy == C_NULL || matrix == C_NULL
        return cholfact_int32(A)
    end

    mknob_out  = ccall((:GetMknob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    dss_handle = ccall((:GetDssHandle, LIB_PATH), Ptr{Void}, (Ptr{Void}, ), mknob_out)
    if dss_handle == NULL
        dss_handle = dss_analyze(structure_proxy)
        ccall((:SetDssHandle, LIB_PATH), Void, (Ptr{Void}, Ptr{Void}), mknob_out, dss_handle)
    end
    dss_factor(dss_handles, matrix)

    # ISSUE: how to return the original result?
end

@doc """ 
y = R \ t

Functionality:
    if fknob == NULL
       return R \ t
   get the input mknob_R from fknob
   if mknob_R->dss_handle == NULL
       return R \ t2
   else
       opt = MKL_DSS_DEFAULTS
       dss_solve!(mknob_R->dss_handle, t2, dy)
"""
function cholmod_factor_inverse_divide(
    y     :: Any,
    R     :: Base.SparseMatrix.CHOLMOD.Factor{Float64},
    t     :: Any,
    fknob :: Ptr{Void}
)
    if fknob == NULL
        y = R \ t
        return
    end

    mknob_in   = ccall((:GetMknob, LIB_PATH), Ptr{Void}, (Ptr{Void}, Cint), fknob, 0)
    dss_handle = ccall((:GetDssHandle, LIB_PATH), Ptr{Void}, (Ptr{Void}, ), mknob_in)
    if dss_handle == NULL
        y = R \ t
        return
    end

    opt = MKL_DSS_DEFAULTS
    dss_solve!(dss_handle, t, y)
end

@doc """
Read a matrix market file. Force it to be symmetric if force_symmetric is true.
Error when symmetric_only is true but the matrix is not symmetric (and not
forced to be symmetric). It calls SpMP library to do the actual work, which is 
faster than a pure Julia version.
"""
function matrix_market_read(
    filename        :: String,
    symmetric_only  :: Bool = false,
    force_symmetric :: Bool = false
)
    # Read into a COO array
    sizes::Vector{Cint} = [0, 0, 0, 0]
    ccall((:load_matrix_market_step1, LIB_PATH), Void, 
        (Ptr{Uint8}, Ptr{Cint}, Bool, Bool), filename, pointer(sizes), force_symmetric, false)

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
        (Ptr{Uint8}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Bool),
        filename, pointer(j), pointer(i), pointer(v), pointer(sizes), true)

    A
end
