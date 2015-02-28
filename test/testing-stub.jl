#TODO: remove this include once the package is ready
include("../src/SparseAccelerator.jl")

using SparseAccelerator
using MatrixMarket

function cg(x, A, b, tol, maxiter)
    r = b - A * x
    rel_err = 1
    p = r
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
    p = z
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

A   = MatrixMarket.mmread("./data/MatrixMarket/BCSSTRUC2/bcsstk14.mtx")

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
#Base.tmerge(Int64, Float64)
#acc_stub(ast[1])
#insert_knobs(ast[1])

@acc x = cg(x, A, b, tol, maxiter)

@acc x = pcg(x, A, b, M, tol, maxiter)

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
