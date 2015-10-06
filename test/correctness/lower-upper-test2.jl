include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)

    trsv_time = 0.
    spmv_time = 0.
    blas1_time = 0.

    L = tril(A)
    U  = spdiagm(1./diag(A))*triu(A)

    for i = 1:2

        B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
        C = tril(A)
        D  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A))*triu(A))
        E = C
        F = D
    end
end

@acc foo()
