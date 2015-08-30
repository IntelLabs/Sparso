include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

# The original pcg_symgs
#include("./pcg-symgs.jl")
function pcg_symgs(x, A, b, tol, maxiter)
    L = tril(A)
    U  = spdiagm(1./diag(A))*triu(A)
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
#        println("iter ", k, ": sum z = ", sum(z), " z=", z)
        
        Base.SparseMatrix.fwdTriSolve!(L, z)
#        println("\tFwdTriSolve! done: sum z = ", sum(z), " z=", z)
        #println("\tL=", L)
          # Could have written as z=L\z if \ is specialized for triangular
          
        #println("\tU before backwadrd=", U)
        Base.SparseMatrix.bwdTriSolve!(U, z)
#        println("\tBwdTriSolve! done: sum z = ", sum(z))#, " z=", z)
        #println("\tU=", U)
          # Could have wrriten as z=U\z if \ is specialized for triangular

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    return x, k, rel_err
end

# This code is what we expect to generate by context-sensitive optimizations.
# It is used for debugging purpose only
function pcg_symgs_with_context_opt(x, A, b, tol, maxiter)
    L = tril(A)
    U  = spdiagm(1./diag(A))*triu(A)
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
    
    __mknobA_8199 = (SparseAccelerator.new_matrix_knob)()
    (SparseAccelerator.set_constant_structured)(__mknobA_8199)
    __fknob_8201 = (SparseAccelerator.new_function_knob)("NewForwardTriangularSolveKnob")
    (SparseAccelerator.add_mknob_to_fknob)(__mknobA_8199,__fknob_8201)
    __fknob_8221 = (SparseAccelerator.new_function_knob)("NewBackwardTriangularSolveKnob")
    (SparseAccelerator.add_mknob_to_fknob)(__mknobA_8199,__fknob_8221)

#println("L=", L)
#println("U=", U)
#println("A=", A)
println("build A1")
tic()
A1 = copy(A)
for i = 1 : size(A1, 1) 
    for j = 1 : size(A1, 2) 
        if i > j
            A1[i, j] = U[j, i]
        end
    end
end
toc()
println("A1 done")
flush(STDOUT)

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
#        println("iter ", k, ": sum z = ", sum(z), " z=", z)
        

        
        z = SparseAccelerator.fwdTriSolve!(L, z, A1,__fknob_8201)
#        println("\tFwdTriSolve! done: sum z = ", sum(z), " z=", z)
        #println("\tL=", L)

        z = SparseAccelerator.bwdTriSolve!(A1,__fknob_8201, z)
#        Base.SparseMatrix.bwdTriSolve!(U, z)
#        println("\tBwdTriSolve! done: sum z = ", sum(z))#, " z=", z)
        #println("\tU=", U)
            
#        z = SparseAccelerator.bwdTriSolve!(A,__fknob_8221, z)
#        println("\tBwdTriSolve! done: sum z = ", sum(z), " z=", z)

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    (SparseAccelerator.delete_function_knob)("DeleteBackwardTriangularSolveKnob",__fknob_8221)
    (SparseAccelerator.delete_function_knob)("DeleteForwardTriangularSolveKnob",__fknob_8201)
    (SparseAccelerator.delete_matrix_knob)(__mknobA_8199)
    
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

