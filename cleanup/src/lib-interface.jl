# This file contains the routines that interfaces with the C library.

# In reordering, we insert some calls to the following 3 functions. So they are executed secretly
# Reorder sparse matrix A and store the result in newA. A itself is not changed.
function CSR_reorder_matrix(A :: SparseMatrixCSC, 
                           newA :: SparseMatrixCSC, 
                           P :: Vector, 
                           Pprime :: Vector, 
                           getPermutation :: Bool, 
                           oneBasedInput :: Bool, 
                           oneBasedOutput :: Bool)
  ccall((:CSR_ReorderMatrix, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool, Bool, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(newA.colptr), pointer(newA.rowval), pointer(newA.nzval),
               pointer(P), pointer(Pprime), getPermutation, oneBasedInput, oneBasedOutput)
end

function CSR_bandwidth(A::SparseMatrixCSC)
   A2 = CreateCSR(A)
   bw = ccall((:CSR_GetBandwidth, LIB_PATH), Cint,
         (Ptr{Void},),
         A2)
   DestroyCSR(A2)
   bw
end

function reorder_vector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVector, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(newV), pointer(P), length(V))
end

function reverse_reorder_vector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVectorWithInversePerm, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(newV), pointer(P), length(V))
end

function CreateCSR(A::SparseMatrixCSC)
  ccall((:CSR_Create, LIB_PATH), Ptr{Void},
       (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint),
       A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval), 1)
end

function DestroyCSR(A::Ptr{Void})
  ccall((:CSR_Destroy, LIB_PATH), Void,
        (Ptr{Void},),
        A)
end

# w = alpha*A*x + beta*y + gamma
function SpMV!(w::Vector, alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector, gamma::Number)
    assert(length(w) == length(y))

    if DEFAULT_LIBRARY == PCL_LIB
        A2 = CreateCSR(A)
        ccall((:CSR_MultiplyWithVector, LIB_PATH), Void,
              (Ptr{Cdouble}, Cdouble, Ptr{Void}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
              pointer(w), alpha, A2, pointer(x), beta, pointer(y), gamma)
        DestroyCSR(A2)
    else
        # use Julia implementation
        w = alpha * A * x + beta * y + gamma
    end
end

# y = A*x
SpMV!(y::Vector, A::SparseMatrixCSC, x::Vector) = SpMV!(y, one(eltype(x)), A, x, zero(eltype(y)), y, zero(eltype(x)))

function PageRank!(w::Vector, alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector, gamma::Number, z::Vector)
    assert(length(w) == length(y))
    assert(length(w) == length(z))

    if DEFAULT_LIBRARY == PCL_LIB
        A2 = CreateCSR(A)
        ccall((:PageRank, LIB_PATH), Void,
              (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
              length(w), pointer(w), alpha, pointer(A.colptr), pointer(A.rowval), pointer(x), beta, pointer(y), gamma, pointer(z))
        DestroyCSR(A2)
    else
        # use Julia implementation
        w = (alpha * A * x + beta * x + gamma).*z
    end
end

# alpha*A*x + beta*y + gamma
function SpMV(alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector, gamma::Number)
  w = Array(Cdouble, length(x))
  SpMV!(w, alpha, A, x, beta, y, gamma)
  w
end

# alpha*A*x + beta*y
SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector, beta::Number, y::AbstractVector) = SpMV(alpha, A, x, beta, y, zero(eltype(x)))

# alpha*A*x + y
SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector, y::AbstractVector) = SpMV(alpha, A, x, one(eltype(y)), y, zero(eltype(x)))

# alpha*A*x
SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector) = SpMV(alpha, A, x, zero(eltype(x)), x, zero(eltype(x)))

# A*x
SpMV(A::SparseMatrixCSC, x::Vector) = SpMV(one(eltype(x)), A, x, zero(eltype(x)), x, zero(eltype(x)))

@doc """
SpMV. In addition, measure the benefit of reordering its input matrix, if
reordering has not done yet.
"""
function SpMV_with_reordering_benefit(A::SparseMatrixCSC, x::Vector, first_reorder_done, beneficial)
    if first_reorder_done
        return SpMV(A, x)
    end
    time1 = time()
    y = SpMV(A, x)
    time2 = time()
    nnz = size(A.nzval, 1)
    rows = A.m
    SpMV_bandwidth = (nnz * 12 + rows * 3 * 8) / (time2 - time1) * 128 / 10 ^ 9
    machine_peak_bandwidth = 60
    if abs(SpMV_bandwidth - machine_peak_bandwidth) > 15
        beneficial[1] = true
    end
    return y
end

function init_conditional_reordering(beneficial, first_reorder_done)
    # Benefical can be just a scalar var. But we will treat it as a 1-element array
    # so that we do not have difficulty in changing it in calling 
    # SpMV_conditional_reordering
    beneficial = Vector{Bool}()
    push!(beneficial, false)
    first_reorder_done = false
end

function WAXPBY!(w::Vector, alpha::Number, x::Vector, beta::Number, y::Vector)
  assert(length(x) == length(y))
  assert(length(x) == length(w))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:waxpby, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
          length(x), pointer(w), alpha, pointer(x), beta, pointer(y))
  else
    w = alpha*x + beta*y
  end
  w
end

function WAXPBY(alpha::Number, x::Vector, beta::Number, y::Vector)
  w = Array(Cdouble, length(x))
  WAXPBY!(w, alpha, x, beta, y)
  w
end

function Dot(x::Vector, y::Vector)
  assert(length(x) == length(y))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:dot, LIB_PATH), Cdouble,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(x), pointer(y))
  else
    dot(x, y)
  end
end

function PointwiseDivide!(w::Vector, x::Vector, y::Vector)
  assert(length(x) == length(y))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:pointwiseDivide, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(w), pointer(x), pointer(y))
  else
    w = x./y
  end
  w
end

function PointwiseDivide(x::Vector, y::Vector)
  w = Array(Cdouble, length(x))
  PointwiseDivide!(w, x, y)
  w
end

function PointwiseMultiply!(w::Vector, x::Vector, y::Vector)
  assert(length(x) == length(y))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:pointwiseMultiply, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(w), pointer(x), pointer(y))
  else
    w = x.*y
  end
  w
end

function PointwiseMultiply(x::Vector, y::Vector)
  w = Array(Cdouble, length(x))
  PointwiseMultiply!(w, x, y)
  w
end

function WAXPB!(w::Vector, alpha::Number, x::Vector, beta::Number)
  assert(length(w) == length(x))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:waxpb, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
          length(x), pointer(w), alpha, pointer(x), beta)
  else
    w = alpha*x + beta
  end
end

function WAXPB(alpha::Number, x::Vector, beta::Number)
  w = Array(Cdouble, length(x))
  WAXPB!(w, alpha, x, beta)
  w
end

end   # end of module

#function Base.A_mul_B!(alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector)
#  SparseAccelerator.SpMV(alpha, A, x, beta, y)
#end

