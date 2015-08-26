include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)
CompilerTools.LivenessAnalysis.set_debug_level(6)

# The original pcg_symgs
include("./pcg-symgs.jl")

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
    
    ___mknobA_8236 = (SparseAccelerator.new_matrix_knob)()
    ___fknob_8238 = (SparseAccelerator.new_function_knob)("NewForwardTriangularSolveKnob")
    (SparseAccelerator.add_mknob_to_fknob)(___mknobA_8236,___fknob_8238)
    ___fknob_8258 = (SparseAccelerator.new_function_knob)("NewBackwardTriangularSolveKnob")
    (SparseAccelerator.add_mknob_to_fknob)(___mknobA_8236,___fknob_8258)
    
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
        
        println("*** Before fwdTriSolve!, A =", A)
        flush(STDOUT)
        
        SparseAccelerator.fwdTriSolve!(A, z, ___fknob_8238)
        SparseAccelerator.bwdTriSolve!(A, z, ___fknob_8258)

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    (SparseAccelerator.delete_function_knob)("DeleteForwardTriangularSolveKnob",___fknob_8238)
    (SparseAccelerator.delete_function_knob)("DeleteBackwardTriangularSolveKnob",___fknob_8258)
    (SparseAccelerator.delete_matrix_knob)(___mknobA_8236) # line 25:
    
    return x, k, rel_err
end

if length(ARGS) == 0
    m = 10
    A = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
else
    A = matrix_market_read(ARGS[1], true, true)
end
m = size(A, 1)
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-10
maxiter = 1000

#println("Original: ")
#@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
#println("\tsum of x=", sum(x))
#println("\tk=", k)
#println("\trel_err=", rel_err)
#flush(STDOUT)

println("\n\nWith manual context-sensitive optimization: ")
x       = zeros(Float64, m)
b       = ones(Float64, m)
x, k, rel_err = pcg_symgs_with_context_opt(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)

