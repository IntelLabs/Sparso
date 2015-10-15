include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(10)

    for i = 1:2
        L, U = SparseAccelerator.ilu(A)
    end
end

@acc foo()
