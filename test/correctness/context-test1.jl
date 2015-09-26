include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

#include("./pcg-symgs.jl")
include("utils.jl")

# TEMPORARY: this should be removed in future once structure analysis can
# recognize that L and U are part of A
# Copied from pcg-symgs.jl to avoid affecting other tests if it were changed.
function pcg_symgs(x, A, b, tol, maxiter)
    # These are the only 3 properties we need to discover automatically.
    # We also require SA_CONST_STRUCTURED for A, but that has already been discovered.
    set_matrix_property(:L, SA_LOWER_OF, :A)
    set_matrix_property(:U, SA_UPPER_OF, :A)
    set_matrix_property(Dict(
        :A => SA_SYMM_STRUCTURED | SA_SYMM_VALUED,
        :U => SA_CONST_VALUED
        )
    )

    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A))*triu(A))
    M = L*U
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1

    z = copy(r)
    Base.SparseMatrix.fwdTriSolve!(L, z)
    Base.SparseMatrix.bwdTriSolve!(U, z)

    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p # Ap = SparseAccelerator.SpMV(A, p) # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        if rel_err < tol 
            break
        end

        z = copy(r)
        Base.SparseMatrix.fwdTriSolve!(L, z)
          # Could have written as z=L\z if \ is specialized for triangular
        Base.SparseMatrix.bwdTriSolve!(U, z)
          # Could have wrriten as z=U\z if \ is specialized for triangular

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    return x, k, rel_err
end

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1], true, true)
    m = size(A, 1)
end
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-10
maxiter = 1000

println("Original: ")
x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tOriginal sum of x=", sum(x))
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

println("\nAccelerated: ")
x = zeros(Float64, m)
@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tAccelerated sum of x=", sum(x))
println("\tAccelerated k=", k)
println("\tAccelerated rel_err=", rel_err)
