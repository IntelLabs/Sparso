include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

function pcg_symgs(x, A, b, tol, maxiter)
    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    blas1_time = 0.

    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)

    spmv_time -= time()
    r = b - A * x
    spmv_time += time()

    blas1_time -= time()
    normr0 = norm(r)
    z = copy(r)
    blas1_time += time()

    rel_err = 1

    trsv_time -= time()
    Base.SparseMatrix.fwdTriSolve!(L, z)
    Base.SparseMatrix.bwdTriSolve!(U, z)
    trsv_time += time()

    blas1_time -= time()
    p = copy(z)
    rz = dot(r, z)
    blas1_time += time()

    k = 1
    while k <= maxiter
        old_rz = rz

        spmv_time -= time()
        Ap = A*p
        spmv_time += time()

        blas1_time -= time()
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        blas1_time += time()

        if rel_err < tol 
            break
        end

        blas1_time -= time()
        z = copy(r)
        blas1_time += time()

        trsv_time -= time()
        L = L \ z
        U = U \ z
        trsv_time += time()

        blas1_time -= time()
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        blas1_time += time()

        k += 1
    end
    total_time = time() - total_time
    println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time")

    return x, k, rel_err
end

if length(ARGS) == 0
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(1000)
else
    # This is for performance testing's purpose
    A = matrix_market_read(ARGS[1], true, true)
end
m       = size(A, 1)
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-7
maxiter = 20000

x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
x = zeros(Float64, m)
x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tOriginal sum of x=", sum(x))
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
x = zeros(Float64, m)
@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tOptimized version: sum of x=", sum(x))
println("\tOptimized version:  k=", k)
println("\tOptimized version: rel_err=", rel_err)