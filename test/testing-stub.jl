#TODO: remove this include once the package is ready
include("../src/OptFramework.jl")
include("../src/SparseAccelerator.jl")
include("../src/sparse-analyze.jl")

using OptFramework
using MatrixMarket

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])

# This function is instrumented with the reordering code we automatically
# generated, so that we can debug the code more easily.
function cg_reordered(x, A, b, tol, maxiter)
    r = b - A * x
    rel_err = 1
    p = copy(r) #NOTE: do not write "p=r"! That would make p and r aliased (the same variable)
    rz = dot(r, r)
    normr0 = sqrt(rz)
    k = 1
    
    __P_51227 = Array(Cint,size(A,2))
    __Pprime_51228 = Array(Cint,size(A,2))
    __A_51229 = SparseMatrixCSC(A.m,A.n,Array(Cint,size(A.colptr,1)),Array(Cint,size(A.rowval,1)),Array(Cdouble,size(A.nzval,1)))
    CSR_ReorderMatrix(A,__A_51229,__P_51227,__Pprime_51228,true)
    A = __A_51229
    __r_51230 = Array(Cdouble,size(r,1))
    reorderVector(r,__r_51230,__P_51227)
    r = __r_51230
    __p_51231 = Array(Cdouble,size(p,1))
    reorderVector(p,__p_51231,__P_51227)
    p = __p_51231
    __x_51232 = Array(Cdouble,size(x,1))
    reorderVector(x,__x_51232,__P_51227)
    x = __x_51232
    
    while k <= maxiter
        old_rz = rz
        Ap = A*p # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rz = dot(r, r)
        rel_err = sqrt(rz)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end
        beta = rz/old_rz
        p = r + beta * p
        k += 1
    end

    __x_51233 = Array(Cdouble,size(x,1))
    reverseReorderVector(x,__x_51233,__P_51227)
    x = __x_51233

    return x, k, rel_err
end

function cg(x, A, b, tol, maxiter)
    r = b - A * x
    rel_err = 1
    p = copy(r) #NOTE: do not write "p=r"! That would make p and r aliased (the same variable)
    rz = dot(r, r)
    normr0 = sqrt(rz)
    k = 1
    while k <= maxiter
        old_rz = rz
        Ap = A*p # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rz = dot(r, r)
        rel_err = sqrt(rz)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end
        beta = rz/old_rz
        p = r + beta * p
        k += 1
    end
    return x, k, rel_err
end


function pcg(x, A, b, M, tol, maxiter)
    r = b - A * x
    z = M \ r
    p = copy(z) # NOTE: do not write "p=z", which make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    while k <= maxiter
        old_rz = rz
        α = old_rz / dot(p, A * p)
        x += α * p
        r -= α * A * p 
        if norm(r) < tol 
            break
        end
        z = M \ r  
        rz = dot(r, z)
        β = rz/old_rz
        p = z + β * p
        k += 1
    end
    return x #, k
end

#A   = MatrixMarket.mmread("./data/MatrixMarket/BCSSTRUC2/bcsstk14.mtx")
# From CSR_test.c in lib/test
A = SparseMatrixCSC{Cdouble, Cint}(
        10,   #m 
        10,   #n
        Cint[1, 3, 7, 10, 14, 17, 21, 25, 27, 28, 29 ],  #colptr
    #int rowPtr[] =    { 0,    2,          6,       9,          13,      16,         20,         24,   26,27, 28 };
        Cint[4, 6, 3, 5, 7, 10, 2, 4, 5, 1, 3, 6, 9, 2, 3, 7, 1, 4, 7, 8, 2, 5, 6, 8, 6, 7, 4, 2], #rowval
    #int colIdx[] =    { 3, 5, 2, 4, 6, 9, 1, 3, 4, 0, 2, 5, 8, 1, 2, 6, 0, 3, 6, 7, 1, 4, 5, 7, 5, 6, 3, 1 };
        
        Cdouble[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]
)

println("******************* A is *************")
println(A)

B = full(A)
println("******************* B is *************")
println(B)


# Julia has only SparseMatrixCSC format so far. But for CG where SpMV is
# important, CSR format is better in performance. However, CSC and CSR
# are the same for symmetric matrices. So for them, we can treat CSC as CSR.
# For a non-symmetric matrix A, to simulate CSR, we
# can transpose A here, and Jongsoo's RCM will treated as a CSR representation,
# even if the matrix type is still "SparseMatrixCSC".
# Note: this actually changes the problem from Ax = b to A.'x = b. 
# A   = A.'

N   = size(A, 1)
M   = speye(N) # Identity
#M   = spdiagm(diag(A)) # Jacobi
#M   = tril(A)*spdiagm(1./diag(A))*triu(A) # Symmetric GS
#M   = A # Perfect

b   = randn(N)
x   = zeros(Float64, N)
tol = 1e-12
maxiter = 2 * N

#println (typeof(x), ", ", typeof(A), ", ",typeof(b), ", ",typeof(M), ", ",typeof(tol), ", ",typeof(maxiter))

#eval(pcg)

#println(names(typeof(methods(pcg, (Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Float64, Int64)))))

#ast = code_lowered(pcg, (Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Float64, Int64))
#println("******************* lowered AST **************")
#println(ast)

#ast = code_typed(pcg, (Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Float64, Int64), optimize=false)
#println("******************* typed AST **************")
#println(ast)

x1, k1, rel_err1 = cg_reordered(x, A, b, tol, maxiter)
println("***** After cg_reordered:")
println("x: ", x1)
println("k: ", k1)
println("rel_err: ", rel_err1)

x   = zeros(Float64, N)
x1, k1, rel_err1 = cg(x, A, b, tol, maxiter)
println("***** After cg:")
println("x: ", x1)
println("k: ", k1)
println("rel_err: ", rel_err1)

#Base.tmerge(Int64, Float64)
#acc_stub(ast[1])
#insert_knobs(ast[1])

@acc x = cg(x, A, b, tol, maxiter)

#@acc x = pcg(x, A, b, M, tol, maxiter)

#println("#Iterations: ", k)
e_v = A*x .- b
abs_err = sqrt(dot(e_v, e_v))
rel_err = abs_err/sqrt(dot(b,b))
#err = sum(abs(A * x .- b))
if (rel_err < 1.0e-6)
	println("Verified")
else
	println("Failed in verification. Error=", rel_err)
end

println("Verifying in Cholfact")
y   = cholfact(A) \ b
e_v = x .- y
abs_err = sqrt(dot(e_v, e_v))
rel_err = abs_err/sqrt(dot(b,b))
#err = sum(abs(x .- y))
println("x .-y = ", rel_err) 
