include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_USE_SPMP, SA_REPLACE_CALLS)

include("utils.jl")

# TEMPORARY: this should be removed in future once structure analysis can
# recognize that L and U are part of A
# Copied from context-test2-ilu0.jl to avoid affecting other tests if it were changed.
function pcg_symgs_ilu0(x, A, b, tol, maxiter)
    Ap = zeros(size(A, 1))

    total_time = time()

    trsv_time = 0.
    spmv_time = 0.
    blas1_time = 0.

    L, U = SparseAccelerator.ilu(A)

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
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    blas1_time += time()

    # For statistics only. We get these values so that L, U and A won't be live out of
    # the loop -- to save reverse-reordering of them.
    nnzL = nnz(L)
    nnzU = nnz(U)
    sizeL1 = size(L, 1) 
    sizeL2 = size(L, 2)
    nnzA = nnz(A)
    sizeA1 = size(A, 1) 
    sizeA2 = size(A, 2)

    k = 1

    while k <= maxiter
        old_rz = rz

        spmv_time -= time()
        Ap = A*p # Ap = SparseAccelerator.SpMV(A, p) # This takes most time. Compiler can reorder A to make faster
        #SparseAccelerator.SpMV!(Ap, A, p)
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
        Base.SparseMatrix.fwdTriSolve!(L, z)
          # Could have written as z=L\z if \ is specialized for triangular
        Base.SparseMatrix.bwdTriSolve!(U, z)
          # Could have wrriten as z=U\z if \ is specialized for triangular
        trsv_time += time()

        blas1_time -= time()
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        blas1_time += time()

        k += 1
    end

    total_time = time() - total_time
    println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnzL   + nnzU  ) + 2.*8*(sizeL1     + sizeL2    ))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnzA   + 8.*(sizeA1     + sizeA2    ))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time")
#   println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time")
 
    return x, k, rel_err
end

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1], true, true)
    m = size(A, 1)
end
x       = zeros(m)
b       = ones(m)
tol     = 1e-7
maxiter = 20000

x, k, rel_err = pcg_symgs_ilu0(x, A, b, tol, maxiter)
x = zeros(m)
println("Original: ")
x, k, rel_err = pcg_symgs_ilu0(x, A, b, tol, maxiter)
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

A2 = copy(A) # workaround that we change A in-place
@acc x, k, rel_err = pcg_symgs_ilu0(x, A2, b, tol, maxiter)
println("\nAccelerated: ")
x = zeros(m)
SparseAccelerator.set_knob_log_level(1)
@acc x, k, rel_err = pcg_symgs_ilu0(x, A, b, tol, maxiter)
println("\tAccelerated k=", k)
println("\tAccelerated rel_err=", rel_err)
