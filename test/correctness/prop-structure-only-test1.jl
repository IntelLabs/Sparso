include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

function foo()
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(10)
    B = spones(A)
    C = speye(20)
    D = B 
    D = A
    E = B
    #D = sparse(A.<3)
end

@acc foo()
