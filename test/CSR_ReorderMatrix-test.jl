using CompilerTools
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using MatrixMarket

include("../Src/SparseAccelerator.jl")
include("./MatrixMarket.jl")
using SparseAccelerator

include("./utils.jl")

println("...... Test 1")

m = 10
A = generate_sparse_matrix(m)
check_symmetry(A)
P = Array(Cint,size(A,2))
Pprime = Array(Cint,size(A,2))
A1 = SparseMatrixCSC(A.m, A.n, Array(Cint, length(A.colptr)), Array(Cint, length(A.rowval)), Array(Cdouble, length(A.nzval)))
SparseAccelerator.CSR_ReorderMatrix(A, A1, P, Pprime, true, true, true)
A = A1

Libc.flush_cstdio()
println("...... Test 1 ends")
flush(STDOUT::IO)

println("...... Test 2") # Extracted from pcg reordered AST

Libc.flush_cstdio()
println("\t Reordering")
flush(STDOUT::IO)

SparseAccelerator.use_lib(SparseAccelerator.PCL_LIB)
A = MatrixMarket.mmread(ARGS[1])
#A = generate_sparse_matrix(m)  # If turning on this, it will pass without errors
check_symmetry(A)
N   = size(A, 1)
b   = ones(Float64, N)
x   = zeros(Float64, N)
inv_d = 1./diag(A)
r = b - A * x
z = inv_d .* r
p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
__P_5162 = Array(Cint,size(A,2))
__Pprime_5163 = Array(Cint,size(A,2))
__A_5164 = SparseMatrixCSC(A.m, A.n, Array(Cint, length(A.colptr)), Array(Cint,length(A.rowval)),Array(Cdouble,length(A.nzval)))
SparseAccelerator.CSR_ReorderMatrix(A,__A_5164,__P_5162,__Pprime_5163,true,true,true)
A = __A_5164
__r_5165 = Array(Cdouble,length(r))
SparseAccelerator.reorderVector(r,__r_5165,__P_5162)
r = __r_5165
__x_5166 = Array(Cdouble,length(x))
SparseAccelerator.reorderVector(x,__x_5166,__P_5162)
x = __x_5166
__p_5167 = Array(Cdouble,length(p))
SparseAccelerator.reorderVector(p,__p_5167,__P_5162)
p = __p_5167

Libc.flush_cstdio()
println("\t Reverse reordering")
flush(STDOUT::IO)

__A_5176 = SparseMatrixCSC(A.m,A.n,Array(Cint,length(A.colptr)),Array(Cint,length(A.rowval)),Array(Cint,length(A.nzval)))
SparseAccelerator.CSR_ReorderMatrix(A,__A_5176,__Pprime_5163,__P_5162,false,true,true)
A = __A_5176

Libc.flush_cstdio()
println("\t\t A reverse reordering done")
flush(STDOUT::IO)


__x_5177 = Array(Cdouble,length(x))
SparseAccelerator.reverseReorderVector(x,__x_5177,__P_5162)::ANY
x = __x_5177
__p_5178 = Array(Cdouble, length(p))
SparseAccelerator.reverseReorderVector(p,__p_5178,__P_5162)
p = __p_5178

Libc.flush_cstdio()
println("...... Test 2 ends")
flush(STDOUT::IO)
