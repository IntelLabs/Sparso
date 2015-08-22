include("../../src/SparseAccelerator.jl")
using SparseAccelerator

# The original pcg_symgs
function pcg_symgs(x, A, b, tol, maxiter)
    L = tril(A)
    U :: SparseMatrixCSC = spdiagm(1./diag(A))*triu(A)
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

        println("\n\nk: ", k)
        println("L: ", L)
        println("U: ", U)
        println("z: ", z)
        flush(STDOUT)

        Base.SparseMatrix.fwdTriSolve!(L, z)
        Base.SparseMatrix.bwdTriSolve!(U, z)

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    return x, k, rel_err
end

# Code that is generated by SparseAccelerator automatically for 
# pcg_symgs(). Used only for manual debugging purpose.
function pcg_symgs_context(x, A, b, tol, maxiter)
    L = tril(A)
    U :: SparseMatrixCSC = spdiagm(1./diag(A))*triu(A)
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

    __mknobL_7876 = SparseAccelerator.new_matrix_knob()
    __mknobU_7878 = SparseAccelerator.new_matrix_knob()
    __fknob_7879 = SparseAccelerator.new_function_knob("NewForwardTriangularSolveKnob")
    SparseAccelerator.add_mknob_to_fknob(__mknobL_7876,__fknob_7879)
    __fknob_7899 = SparseAccelerator.new_function_knob("NewBackwardTriangularSolveKnob")
    SparseAccelerator.add_mknob_to_fknob(__mknobU_7878,__fknob_7899)

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

        println("\n\nk: ", k)
        println("L: ", L)
        println("U: ", U)
        println("z: ", z)
        flush(STDOUT)

        SparseAccelerator.fwdTriSolve!(L, z, __fknob_7879)
        SparseAccelerator.bwdTriSolve!(U, z, __fknob_7899)

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    
    SparseAccelerator.delete_function_knob("DeleteForwardTriangularSolveKnob",__fknob_7879)
    SparseAccelerator.delete_function_knob("DeleteBackwardTriangularSolveKnob",__fknob_7899)
    SparseAccelerator.delete_matrix_knob(__mknobL_7876)
    SparseAccelerator.delete_matrix_knob(__mknobU_7878)
    
    return x, k, rel_err
end

A = matrix_market_read(ARGS[1])
m = size(A, 1)
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

println("\n\nWith context: ")
x       = zeros(Float64, m)
b       = ones(Float64, m)
x, k, rel_err = pcg_symgs_context(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)

