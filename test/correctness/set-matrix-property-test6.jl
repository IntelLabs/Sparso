include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    set_matrix_property(Dict(
        :A => SA_STRUCTURE_ONLY, 
        )
    )

    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    B = generate_symmetric_nonzerodiagonal_sparse_matrix(m)

    for i = 1:2
        set_matrix_property(Dict(
            :C => SA_STRUCTURE_ONLY, 
            )
        )

        C = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
        D = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    end
end

@acc foo()
