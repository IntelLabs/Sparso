include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    set_matrix_property(Dict(
        :S => SA_SYMM_VALUED, 
        )
    )

    m = 10
    S = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    NS = generate_symmetric_nonzerodiagonal_sparse_matrix(m)

    for i = 1:2
        B = S + S + S # symm
        C = S + S + NS
        D = S - NS
        E = B + B # symm
        F = E - S # symm
        G = S + B - F # symm
        H = NS - B + G
        G = G + G - G
    end
end

@acc foo()
