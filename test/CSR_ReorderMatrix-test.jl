include("../src/SparseAccelerator.jl")
using SparseAccelerator

include("./utils.jl")

m = 10
A = generate_sparse_matrix(m)
check_symmetry(A)

P = Array(Cint,size(A,2))
Pprime = Array(Cint,size(A,2))
A1 = SparseMatrixCSC(A.m, A.n, Array(Cint, length(A.colptr)), Array(Cint, length(A.rowval)), Array(Cdouble, length(A.nzval)))

SparseAccelerator.CSR_ReorderMatrix(A, A1, P, Pprime, true, true, true)
A = A1
