include("../src/SparseAccelerator.jl")
using SparseAccelerator

m = 10
A = sprand(m, m, 0.1)
for i = 1:m for j = 1:m A[i, j] = 0.0 end end
A[6 ,  1]  =  0.662534
A[10,  1]  =  0.123341
A[4 ,  4]  =  0.784948
A[8 ,  4]  =  0.993791
A[10,  5]  =  0.283148
A[1 ,  6]  =  0.662534
A[4 ,  8]  =  0.993791
A[8 ,  8]  =  0.969205
A[1 , 10]  =  0.123341
A[5 , 10]  =  0.283148

P = Array(Cint,size(A,2))
Pprime = Array(Cint,size(A,2))
A1 = SparseMatrixCSC(A.m, A.n, Array(Cint, length(A.colptr)), Array(Cint, length(A.rowval)), Array(Cdouble, length(A.nzval)))

println("**** Checking symmetry")
for i = 1:m 
    for j = 1:m 
        if A[i, j] != A[j, i]
            println("Matrix is asymmetric!")
            println("A[", i, ",", j, "] != A[", j, ",", i, "]")
        end
    end
end
println("Done checking")
flush(STDOUT::IO)

SparseAccelerator.CSR_ReorderMatrix(A, A1, P, Pprime, true, true, true)
A = A1
