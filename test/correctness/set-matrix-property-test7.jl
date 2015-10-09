include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS)

include("utils.jl")

function foo()
    set_matrix_property(:A, SA_TRANSPOSE_OF, :B)

    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)

    for i = 1:m
        set_matrix_property(:B, SA_TRANSPOSE_OF, :C)

        C = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
        D = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    end
end

@acc foo()
