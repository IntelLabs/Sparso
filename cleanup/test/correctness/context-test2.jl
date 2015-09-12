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
    #r = b - A * x
    r = b - SparseAccelerator.SpMV(A,x)
    spmv_time += time()

    blas1_time -= time()
    #normr0 = norm(r)
    normr0 = SparseAccelerator.norm(r)
    z = copy(r)
    blas1_time += time()

    rel_err = 1

    trsv_time -= time()
    Base.SparseMatrix.fwdTriSolve!(L, z)
    Base.SparseMatrix.bwdTriSolve!(U, z)
    trsv_time += time()

    blas1_time -= time()
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    #rz = dot(r, z)
    rz = SparseAccelerator.dot(r,z)
    blas1_time += time()

    k = 1
    while k <= maxiter
        old_rz = rz

        spmv_time -= time()
        #Ap = A*p
        Ap = SparseAccelerator.SpMV(A,p)
        spmv_time += time()

        blas1_time -= time()
        #alpha = old_rz / dot(p, Ap)
        alpha = old_rz / SparseAccelerator.dot(p, Ap)
        #x += alpha * p
        SparseAccelerator.WAXPBY!(x,1,x,alpha,p)
        #r -= alpha * Ap
        SparseAccelerator.WAXPBY!(r,1,r,-alpha,Ap)
        #rel_err = norm(r)/normr0
        rel_err = SparseAccelerator.norm(r)/normr0
        blas1_time += time()

        if rel_err < tol 
            break
        end

        blas1_time -= time()
        #z = copy(r)
        SparseAccelerator.copy!(z, r)
        blas1_time += time()

        trsv_time -= time()
        Base.SparseMatrix.fwdTriSolve!(L, z)
        Base.SparseMatrix.bwdTriSolve!(U, z)
        trsv_time += time()

        blas1_time -= time()
        #rz = dot(r, z)
        rz = SparseAccelerator.dot(r,z)
        beta = rz/old_rz
        #p = z + beta * p
        SparseAccelerator.WAXPBY!(p,1,z,beta,p)
        blas1_time += time()

        k += 1
    end
    total_time = time() - total_time
    println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time")

    return x, k, rel_err
end

# This code is what we expect to generate by context-sensitive optimizations except for reordering.
# It is used for debugging purpose only
function pcg_symgs_with_context_opt_without_reordering(x, A, b, tol, maxiter)
    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    dot_time = 0.
    blas1_time = 0.

    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)

    spmv_time -= time()
    #r = b - A * x
    r = b - SparseAccelerator.SpMV(A,x)
    spmv_time += time()

    blas1_time -= time()
    #normr0 = norm(r)
    normr0 = SparseAccelerator.norm(r)
    z = copy(r)
    blas1_time += time()

    rel_err = 1

    trsv_time -= time()
    Base.SparseMatrix.fwdTriSolve!(L, z)
    Base.SparseMatrix.bwdTriSolve!(U, z)
    trsv_time += time()

    blas1_time -= time()
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    #rz = dot(r, z)
    rz = SparseAccelerator.dot(r,z)
    blas1_time += time()

    k = 1

    __mknobA = (SparseAccelerator.new_matrix_knob)(A, true, true, true, true, false, false)
    __mknobL = (SparseAccelerator.new_matrix_knob)(L, true, true, false, false, false, false) # matrix knob for L
    __mknobU = (SparseAccelerator.new_matrix_knob)(U, true, true, false, false, false, false) # matrix knob for L

    (SparseAccelerator.set_derivative)(__mknobL, SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)
    (SparseAccelerator.set_derivative)(__mknobU, SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)

    __fknob_8201 = (SparseAccelerator.new_function_knob)()
    (SparseAccelerator.add_mknob_to_fknob)(__mknobL,__fknob_8201)
    __fknob_8221 = (SparseAccelerator.new_function_knob)()
    (SparseAccelerator.add_mknob_to_fknob)(__mknobU,__fknob_8221)

    while k <= maxiter
        old_rz = rz

        spmv_time -= time()
        #Ap = A*p
        Ap = SparseAccelerator.SpMV(A,p)
        spmv_time += time()

        blas1_time -= time()
        #alpha = old_rz / dot(p, Ap)
        alpha = old_rz / SparseAccelerator.dot(p, Ap)
        #x += alpha * p
        SparseAccelerator.WAXPBY!(x,1,x,alpha,p)
        #r -= alpha * Ap
        SparseAccelerator.WAXPBY!(r,1,r,-alpha,Ap)
        #rel_err = norm(r)/normr0
        rel_err = SparseAccelerator.norm(r)/normr0
        blas1_time += time()

        if rel_err < tol 
            break
        end

        blas1_time -= time()
        #z = copy(r)
        SparseAccelerator.copy!(z, r)
        blas1_time += time()

        trsv_time -= time()
        #Base.SparseMatrix.fwdTriSolve!(L, z)
        SparseAccelerator.fwdTriSolve!(L, z, __fknob_8201)
        #Base.SparseMatrix.bwdTriSolve!(U, z)
        SparseAccelerator.bwdTriSolve!(U, z, __fknob_8221)
        trsv_time += time()

        blas1_time -= time()
        #rz = dot(r, z)
        rz = SparseAccelerator.dot(r,z)
        beta = rz/old_rz
        #p = z + beta * p
        SparseAccelerator.WAXPBY!(p,1,z,beta,p)
        blas1_time += time()

        k += 1
    end

    (SparseAccelerator.delete_function_knob)(__fknob_8221)
    (SparseAccelerator.delete_function_knob)(__fknob_8201)
    (SparseAccelerator.delete_matrix_knob)(__mknobL)
    
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
    reorder_time = 0.

    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)

    spmv_time -= time()
    #r = b - A * x
    r = b - SparseAccelerator.SpMV(A,x)
    spmv_time += time()

    blas1_time -= time()
    #normr0 = norm(r)
    normr0 = SparseAccelerator.norm(r)
    z = copy(r)
    blas1_time += time()

    rel_err = 1

    trsv_time -= time()
    Base.SparseMatrix.fwdTriSolve!(L, z)
    Base.SparseMatrix.bwdTriSolve!(U, z)
    trsv_time += time()

    blas1_time -= time()
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    #rz = dot(r, z)
    rz = SparseAccelerator.dot(r,z)
    blas1_time += time()

    k = 1

    __mknobA = (SparseAccelerator.new_matrix_knob)(A, true, true, true, true, false, false)
    __mknobL = (SparseAccelerator.new_matrix_knob)(L, true, true, false, false, false, false) # matrix knob for L
    __mknobU = (SparseAccelerator.new_matrix_knob)(U, true, true, false, false, false, false) # matrix knob for L

    (SparseAccelerator.set_derivative)(__mknobL, SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)
    (SparseAccelerator.set_derivative)(__mknobU, SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)

    __fknob_8201 = (SparseAccelerator.new_function_knob)()
    (SparseAccelerator.add_mknob_to_fknob)(__mknobL,__fknob_8201)
    __fknob_8221 = (SparseAccelerator.new_function_knob)()
    (SparseAccelerator.add_mknob_to_fknob)(__mknobU,__fknob_8221)

    # Set up reordering:
    # Let SparseAccelerator.fwdTriSolve! decides what permutation/inverse
    # permutation vector should be, and reorders its inputs andoutputs accordingly.
    # Other functions just respect the decision, makes no decision, nor does any
    # reordering.
    (SparseAccelerator.set_reordering_decision_maker)(__fknob_8201)

    perm                         = C_NULL
    inv_perm                     = C_NULL
    reordering_on_back_edge_done = false
    initial_reordering_done      = false
    while k <= maxiter
        if initial_reordering_done && !reordering_on_back_edge_done
            # This is the first time the execution reaches here (along the 
            # backedge of the loop) after the initial reordering was done.
            # Reorder all the arrays that are affected by the other
            # already-reordered arrays.
            assert(perm != C_NULL && inv_perm != C_NULL)

            reorder_time -= time()

            new_A = copy(A)
            SparseAccelerator.reorder_matrix(A, new_A, perm, inv_perm)
            A = new_A

            reorder_time += time()

            reordering_on_back_edge_done = true
        end
    
    
        old_rz = rz

        spmv_time -= time()
        #Ap = A*p
        Ap = SparseAccelerator.SpMV(A,p)
        spmv_time += time()

        blas1_time -= time()
        #alpha = old_rz / dot(p, Ap)
        alpha = old_rz / SparseAccelerator.dot(p, Ap)
        #x += alpha * p
        SparseAccelerator.WAXPBY!(x,1,x,alpha,p)
        #r -= alpha * Ap
        SparseAccelerator.WAXPBY!(r,1,r,-alpha,Ap)
        #rel_err = norm(r)/normr0
        rel_err = SparseAccelerator.norm(r)/normr0
        blas1_time += time()

        if rel_err < tol 
            break
        end

        blas1_time -= time()
        #z = copy(r)
        SparseAccelerator.copy!(z, r)
        blas1_time += time()

        trsv_time -= time()
        #Base.SparseMatrix.fwdTriSolve!(L, z)
        SparseAccelerator.fwdTriSolve!(L, z, __fknob_8201)
        trsv_time += time()

        if !initial_reordering_done
            reorder_time -= time()

            perm, inv_perm = SparseAccelerator.get_reordering_vectors(__fknob_8201)
            if perm != C_NULL
                assert(inv_perm != C_NULL)

                # This is the first time fwdTriSolve! has decided to do reordering.
                # Its own inputs and outputs (L and z) have already been reordered.
                # We need to reorder other arrays that affected by L and z in the
                # subsequent execution.
                new_U = copy(U)
                SparseAccelerator.reorder_matrix(U, new_U, perm, inv_perm)
                U = new_U
                
                new_r = copy(r)
                SparseAccelerator.reorder_vector(r, new_r, perm)
                r = new_r
                
                new_p = copy(p)
                SparseAccelerator.reorder_vector(p, new_p, perm)
                p = new_p

                initial_reordering_done = true
                # For any static program point from here to the end of the loop,
                # any array that needs to be reordered has been reordered.
            end

            reorder_time += time()
        end

        trsv_time -= time()
        #Base.SparseMatrix.bwdTriSolve!(U, z)
        SparseAccelerator.bwdTriSolve!(U, z, __fknob_8221)
        trsv_time += time()

        blas1_time -= time()
        #rz = dot(r, z)
        rz = SparseAccelerator.dot(r,z)
        beta = rz/old_rz
        #p = z + beta * p
        SparseAccelerator.WAXPBY!(p,1,z,beta,p)
        blas1_time += time()

        k += 1
    end

    if reordering_on_back_edge_done
        # Reverse reorder live arrays that have been reordered.
        # Only x will live out here, and it is reordered only if reordering_on_back_edge_done
        reorder_time -= time()
        new_x = copy(x)
        SparseAccelerator.reverse_reorder_vector(x, new_x, perm)
        x = new_x
        reorder_time += time()
    end
    
    (SparseAccelerator.delete_function_knob)(__fknob_8221)
    (SparseAccelerator.delete_function_knob)(__fknob_8201)
    (SparseAccelerator.delete_matrix_knob)(__mknobL)
    
    total_time = time() - total_time
    println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time reorder_time = $reorder_time")

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
tol     = 1e-5
maxiter = 1000

println("Original: ")
x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)
flush(STDOUT)

println("\nWith manual context-sensitive optimization without reordering: ")
x       = zeros(Float64, m)
b       = ones(Float64, m)
x, k, rel_err = pcg_symgs_with_context_opt_without_reordering(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)

println("\nWith manual context-sensitive optimization: ")
x       = zeros(Float64, m)
b       = ones(Float64, m)
x, k, rel_err = pcg_symgs_with_context_opt(x, A, b, tol, maxiter)
println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)

