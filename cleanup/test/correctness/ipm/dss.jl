const LIB_CSR_PATH = "../../../../lib/libcsr.so"

const MKL_DSS_DEFAULTS = 0
const MKL_DSS_NON_SYMMETRIC = 536871104
const MKL_DSS_SUCCESS = 0
const MKL_DSS_AUTO_ORDER = 268435520
const MKL_DSS_POSITIVE_DEFINITE = 134217792

function dss_analyze(A::SparseMatrixCSC)
  handle = Int[0]
  opt = MKL_DSS_DEFAULTS
  error = ccall((:dss_create, LIB_CSR_PATH), Cint,
        (Ptr{Void}, Ptr{Cint}),
        handle, &opt)
  if error != MKL_DSS_SUCCESS
    println("dss_create returned error code $error")
  end

  opt = MKL_DSS_NON_SYMMETRIC
  m = size(A, 1)
  nnz = A.colptr[m + 1] - 1
  error = ccall((:dss_define_structure, LIB_CSR_PATH), Cint,
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        handle, &opt, A.colptr, &m, &m, A.rowval, &nnz)
  if error != MKL_DSS_SUCCESS
    println("dss_define_structure returned error code $error")
  end

  opt = MKL_DSS_AUTO_ORDER
  error = ccall((:dss_reorder, LIB_CSR_PATH), Cint,
                (Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
                handle, &opt, C_NULL)
  if error != MKL_DSS_SUCCESS
    println("dss_reorder returned error code $error")
  end

  handle
end

function dss_factor(handle, A::SparseMatrixCSC)
  opt = MKL_DSS_POSITIVE_DEFINITE
  error = ccall((:dss_factor_real, LIB_CSR_PATH), Cint,
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
  error = ccall((:dss_solve_real, LIB_CSR_PATH), Cint,
                (Ptr{Void}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
                handle, &opt, rhs, &nrhs, sol)
end
