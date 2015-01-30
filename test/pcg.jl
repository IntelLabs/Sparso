using SparseAccelerator
using MatrixMarket

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
    return x, k
end

A   = MatrixMarket.mmread("./data/MatrixMarket/BCSSTRUC2/bcsstk14.mtx")
N   = size(A, 1)
M   = speye(N) # Identity
#M   = spdiagm(diag(A)) # Jacobi
#M   = tril(A)*spdiagm(1./diag(A))*triu(A) # Symmetric GS
#M   = A # Perfect

b   = randn(N)
x   = zeros(Float64, N)
tol = 1e-12
maxiter = 2 * N
@acc x, k = pcg(x, A, b, M, tol, maxiter)

println("#Iterations: ", k)
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
