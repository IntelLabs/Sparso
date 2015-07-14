using CompilerTools
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis
using MatrixMarket

include("../Src/SparseAccelerator.jl")
include("./MatrixMarket.jl")
using SparseAccelerator

SparseAccelerator.use_lib(SparseAccelerator.PCL_LIB)
A = MatrixMarket.mmread(ARGS[1])

Libc.flush_cstdio()
println("\t Reordering")
flush(STDOUT::IO)

__P_5162 = Array(Cint,size(A,2))
__Pprime_5163 = Array(Cint,size(A,2))
__A_5164 = SparseMatrixCSC(A.m, A.n, Array(Cint, length(A.colptr)), Array(Cint,length(A.rowval)),Array(Cdouble,length(A.nzval)))
SparseAccelerator.CSR_ReorderMatrix(A,__A_5164,__P_5162,__Pprime_5163,true,true,true)
A = __A_5164

Libc.flush_cstdio()
println("\t Reverse reordering")
flush(STDOUT::IO)

__A_5176 = SparseMatrixCSC(A.m,A.n,Array(Cint,length(A.colptr)),Array(Cint,length(A.rowval)),Array(Cint,length(A.nzval)))
SparseAccelerator.CSR_ReorderMatrix(A,__A_5176,__Pprime_5163,__P_5162,false,true,true)
A = __A_5176
