include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(10)
    At = A'
    Att = At'
    B = At
    C = Att
    C = B
end

@acc foo()
