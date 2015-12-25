# $ julia knob-test1.jl ../matrices/bcsstk14.mtx 
# assymmetric (2, 8) exists but (8, 2) doesn't
# julia: SymGS.cpp:305: bool SpMP::getSymmetricNnzPattern(const SpMP::CSR *, int **, int **, int **, int **): Assertion `sym.isSymmetric(false, true)' failed.
# Aborted (core dumped)

include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

function pcg_symgs(x, A, b)
    L, U = SparseAccelerator.ilu(A)
    r = b - SparseAccelerator.SpMV(1,A,x)
    z = copy(r)
    SparseAccelerator.fwdTriSolve!(z,L,r, C_NULL)
end

A = matrix_market_read(ARGS[1], true, true)
m = size(A, 1)
b = ones(Float64, m)
x = zeros(Float64, m)
pcg_symgs(x, A, b)