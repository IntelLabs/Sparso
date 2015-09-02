include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

# The original pcg_symgs
#include("./pcg-symgs.jl")
function pcg_symgs(x, A, b, tol, maxiter)
    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    dot_time = 0.
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
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
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
        Base.SparseMatrix.fwdTriSolve!(L, z)
        Base.SparseMatrix.bwdTriSolve!(U, z)
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

# This code is what we expect to generate by context-sensitive optimizations.
# It is used for debugging purpose only
function pcg_symgs_with_context_opt(x, A, b, tol, maxiter)
    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    dot_time = 0.
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
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    blas1_time += time()

    k = 1

    __mknobA = (SparseAccelerator.new_matrix_knob)(A)
    __mknobL = (SparseAccelerator.new_matrix_knob)(L) # matrix knob for L
    __mknobU = (SparseAccelerator.new_matrix_knob)(U) # matrix knob for L

    (SparseAccelerator.set_derivative)(__mknobL, SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)
    (SparseAccelerator.set_derivative)(__mknobU, SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)

    (SparseAccelerator.set_constant_valued)(__mknobA)
    (SparseAccelerator.set_value_symmetric)(__mknobA)
    (SparseAccelerator.set_constant_valued)(__mknobL)
    (SparseAccelerator.set_constant_valued)(__mknobU)
    __fknob_8201 = (SparseAccelerator.new_function_knob)("NewForwardTriangularSolveKnob")
    (SparseAccelerator.add_mknob_to_fknob)(__mknobL,__fknob_8201)
    __fknob_8221 = (SparseAccelerator.new_function_knob)("NewBackwardTriangularSolveKnob")
    (SparseAccelerator.add_mknob_to_fknob)(__mknobU,__fknob_8221)

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
        SparseAccelerator.fwdTriSolve!(L, z, __fknob_8201)
        SparseAccelerator.bwdTriSolve!(U, z, __fknob_8221)
        trsv_time += time()

        blas1_time -= time()
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        blas1_time += time()

        k += 1
    end
    (SparseAccelerator.delete_function_knob)("DeleteBackwardTriangularSolveKnob",__fknob_8221)
    (SparseAccelerator.delete_function_knob)("DeleteForwardTriangularSolveKnob",__fknob_8201)
    (SparseAccelerator.delete_matrix_knob)(__mknobL)
    
    total_time = time() - total_time
    println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time")

    return x, k, rel_err
end

if length(ARGS) == 0
    m = 1000
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
#    A = matrix_market_read("tiny-diag.mtx", true, true)
else
    A = matrix_market_read(ARGS[1], true, true)
end
m       = size(A, 1)
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-10
maxiter = 1000

println("Original: ")
x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)
flush(STDOUT)

println("\n\nWith manual context-sensitive optimization: ")
x       = zeros(Float64, m)
b       = ones(Float64, m)
x, k, rel_err = pcg_symgs_with_context_opt(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)

