# Note: tested with julia-b77587362d (0.5.0 dev) nightly build. Julia-0.4 has
# issues with "\", and cannot run it through.
 
include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS, SA_REORDER)

function pcg_symgs(x, A, b, tolerance, maxiter)
    set_matrix_property(:L, SA_LOWER_OF, :A)
    set_matrix_property(:U, SA_UPPER_OF, :A)
    set_matrix_property(Dict(
        :A => SA_SYMM_STRUCTURED | SA_SYMM_VALUED,
        )
    )

    total_time = time()

    # TODO: replace this ILU0 with ICHOL
    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)
    r = b - A * x
    normr0 = norm(r)
    z = L \ r    
    z = U \ z    
    p = copy(z)
    rz = dot(r, z)

    rel_err = 1
    k = 1
    while true
        old_rz = rz
        Ap = A*p
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        if rel_err < tolerance 
            break
        end
        z = L \ r
        z = U \ z
        rz = dot(r, z)
        p = z + (rz/old_rz) * p
        k += 1
    end

    total_time = time() - total_time
    println("total = $(total_time)s")

    return x, k, rel_err
end

if length(ARGS) == 0
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(1000)
else
    # This is for performance testing's purpose
    A = matrix_market_read(ARGS[1], true, true)
end
m         = size(A, 1)
b         = ones(Float64, m)
originalA = copy(A)
tolerance = 1e-7
maxiter   = 20000

# Run each code version twice. Measure the time for the second time.
x = zeros(Float64, m)
x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
x = zeros(Float64, m)
x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
println("\tOriginal sum of x=", sum(x))
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

x = zeros(Float64, m)
@acc x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
x = zeros(Float64, m)
b = ones(Float64, m)
A = originalA
@acc x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
println("\tOptimized version: sum of x=", sum(x))
println("\tOptimized version:  k=", k)
println("\tOptimized version: rel_err=", rel_err)