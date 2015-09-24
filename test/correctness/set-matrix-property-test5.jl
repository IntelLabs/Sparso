include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    set_matrix_property(:A, SA_LOWER_OF, :B)

    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)

    for i = 1:2
        B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
        m = size(A, 1)
        B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    end

    for i = 1:2
        set_matrix_property(:B, SA_UPPER_OF, :C)
        B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
        m = size(A, 1)
        C = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    end
end

@acc foo()
