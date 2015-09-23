include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    set_matrix_property(Dict(
        :A => SA_CONST_VALUED | SA_SYMM_VALUED | SA_SYMM_STRUCTURED, 
        :B => SA_CONST_VALUED | SA_CONST_STRUCTURED | SA_SYMM_VALUED
        )
    )

    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)


    for i = 1:2
        B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
        m = size(A, 1)
        B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    end
end

@acc foo()