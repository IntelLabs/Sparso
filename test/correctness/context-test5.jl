#=
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

include("../../src/Sparso.jl")
include("../../src/simple-show.jl")
using Sparso

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)

include("utils.jl")

# TEMPORARY: this should be removed in future once structure analysis can
# recognize that L and U are part of A
# Copied from context-test2-ilu0.jl to avoid affecting other tests if it were changed.
function pcg_symgs_ilu0(x, A, b, tol, maxiter)
    # These are the only 3 properties we need to discover automatically.
    # We also require SA_CONST_STRUCTURED for A, but that has already been discovered.
    set_matrix_property(Dict(
        :A => SA_SYMM_VALUED
        )
    )

    Ap = zeros(size(A, 1))

    total_time = time()

    trsv_time = 0.
    spmv_time = 0.
    blas1_time = 0.

    L, U = Sparso.ilu(A)

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

    bytesTrsvPerIter = 12.*(nnz(L)+nnz(U)) + 2.*4*size(L,1) + 2.*8*(size(L,1) + size(L,2))
    bytesSpMVPerIter = 12.*nnz(A) + 4*size(A,1) + 8*(size(A,1) + size(A,2))

    k = 1

    # Pre-allocation of space. Should get rid of once Linxiang's object removal works
    # TODO: remove it.
    Ap = Array(Cdouble, length(p))
    while k <= maxiter
        set_matrix_property(Dict(
            :r  => SA_HAS_DEDICATED_MEMORY,
            :x  => SA_HAS_DEDICATED_MEMORY,
            :p  => SA_HAS_DEDICATED_MEMORY,
          )
        )
    
        old_rz = rz

        spmv_time -= time()
        #Ap = A*p # Ap = Sparso.SpMV(A, p) # This takes most time. Compiler can reorder A to make faster
        #Sparso.SpMV!(Ap, A, p)
        # TODO: remove A_mul_B!. Use A*p once object removal works
        A_mul_B!(Ap, A, p)
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

    bytesTrsv = bytesTrsvPerIter*k
    bytesSpMV = bytesSpMVPerIter*(k + 1)

    println("total = $(total_time)s ($((bytesTrsv + bytesSpMV)/total_time/1e9) gbps) trsv_time = $(trsv_time)s ($(bytesTrsv/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($(bytesSpMV/spmv_time/1e9) gbps) blas1_time = $blas1_time")
 
    return x, k, rel_err
end

# This code is what we expect to generate by context-sensitive optimizations except for reordering.
# It is used for debugging purpose only
function pcg_symgs_ilu0_with_context_opt_without_reordering(x, A, b, tol, maxiter, do_print)
    Ap = zeros(size(A, 1))

    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    dot_time = 0.
    blas1_time = 0.

    L, U = Sparso.ilu(A)

    spmv_time -= time()
    #r = b - A * x
    r = b - Sparso.SpMV(A,x)
    spmv_time += time()

    blas1_time -= time()
    #normr0 = norm(r)
    normr0 = Sparso.norm(r)
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
    rz = Sparso.dot(r,z)
    blas1_time += time()

    k = 1

    __mknobA = (Sparso.new_matrix_knob)(:A, true, true, true, true, false, false)
    __mknobL = (Sparso.new_matrix_knob)(:L, true, true, false, false, false, false) # matrix knob for L
    __mknobU = (Sparso.new_matrix_knob)(:U, true, true, false, false, false, false) # matrix knob for L

    (Sparso.set_derivative)(__mknobL, Sparso.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)
    (Sparso.set_derivative)(__mknobU, Sparso.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)

    fknob_spmv = (Sparso.new_function_knob)()
    (Sparso.add_mknob_to_fknob)(__mknobA,fknob_spmv)
    __fknob_8201 = (Sparso.new_function_knob)()
    (Sparso.add_mknob_to_fknob)(__mknobL,__fknob_8201)
    __fknob_8221 = (Sparso.new_function_knob)()
    (Sparso.add_mknob_to_fknob)(__mknobU,__fknob_8221)

    while k <= maxiter
        old_rz = rz

        spmv_time -= time()
        #Ap = A*p
        Sparso.SpMV!(Ap, A,p,fknob_spmv)
        spmv_time += time()

        blas1_time -= time()
        #alpha = old_rz / dot(p, Ap)
        alpha = old_rz / Sparso.dot(p, Ap)
        #x += alpha * p
        Sparso.WAXPBY!(x,1,x,alpha,p)
        #r -= alpha * Ap
        Sparso.WAXPBY!(r,1,r,-alpha,Ap)
        #rel_err = norm(r)/normr0
        rel_err = Sparso.norm(r)/normr0
        blas1_time += time()

        if rel_err < tol 
            break
        end

        blas1_time -= time()
        #z = copy(r)
        Sparso.copy!(z, r)
        blas1_time += time()

        trsv_time -= time()
        #Base.SparseMatrix.fwdTriSolve!(L, z)
        Sparso.fwdTriSolve!(L, z, __fknob_8201)
        #Base.SparseMatrix.bwdTriSolve!(U, z)
        Sparso.bwdTriSolve!(U, z, __fknob_8221)
        trsv_time += time()

        blas1_time -= time()
        #rz = dot(r, z)
        rz = Sparso.dot(r,z)
        beta = rz/old_rz
        #p = z + beta * p
        Sparso.WAXPBY!(p,1,z,beta,p)
        blas1_time += time()

        k += 1
    end

    (Sparso.delete_function_knob)(__fknob_8221)
    (Sparso.delete_function_knob)(__fknob_8201)
    (Sparso.delete_matrix_knob)(__mknobL)
    
    if do_print
      total_time = time() - total_time
      println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time")
    end

    return x, k, rel_err
end

# This code is what we expect to generate by context-sensitive optimizations.
# It is used for debugging purpose only
function pcg_symgs_ilu0_with_context_opt(x, A, b, tol, maxiter, do_print)
    A = copy(A)
    Ap = zeros(size(A, 1))

    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    dot_time = 0.
    blas1_time = 0.
    reorder_time = 0.

    L, U = Sparso.ilu(A)

    spmv_time -= time()
    #r = b - A * x
    r = b - Sparso.SpMV(A,x)
    spmv_time += time()

    blas1_time -= time()
    #normr0 = norm(r)
    normr0 = Sparso.norm(r)
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
    rz = Sparso.dot(r,z)
    blas1_time += time()

    k = 1

    __mknobA = (Sparso.new_matrix_knob)(:A, true, true, true, true, false, false)
    __mknobL = (Sparso.new_matrix_knob)(:L, true, true, false, false, false, false) # matrix knob for L
    __mknobU = (Sparso.new_matrix_knob)(:U, true, true, false, false, false, false) # matrix knob for L

    (Sparso.set_derivative)(__mknobL, Sparso.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)
    (Sparso.set_derivative)(__mknobU, Sparso.DERIVATIVE_TYPE_SYMMETRIC, __mknobA)

    fknob_spmv = (Sparso.new_function_knob)()
    (Sparso.add_mknob_to_fknob)(__mknobA,fknob_spmv)
    __fknob_8201 = (Sparso.new_function_knob)()
    (Sparso.add_mknob_to_fknob)(__mknobL,__fknob_8201)
    __fknob_8221 = (Sparso.new_function_knob)()
    (Sparso.add_mknob_to_fknob)(__mknobU,__fknob_8221)

    # Set up reordering:
    # Let Sparso.fwdTriSolve! decides what permutation/inverse
    # permutation vector should be, and reorders its inputs andoutputs accordingly.
    # Other functions just respect the decision, makes no decision, nor does any
    # reordering.
    (Sparso.set_reordering_decision_maker)(__fknob_8201)

    # Reordering_status is a vector. See the comments in lib-interface.jl:
    # reordering() for the meaning of the elements
    reordering_status = [false, C_NULL, C_NULL, C_NULL, C_NULL, reorder_time]
    while k <= maxiter
        old_rz = rz

        spmv_time -= time()
        #Ap = A*p
        Sparso.SpMV!(Ap,A,p,fknob_spmv)
        spmv_time += time()

        blas1_time -= time()
        #alpha = old_rz / dot(p, Ap)
        alpha = old_rz / Sparso.dot(p, Ap)
        #x += alpha * p
        Sparso.WAXPBY!(x,1,x,alpha,p)
        #r -= alpha * Ap
        Sparso.WAXPBY!(r,1,r,-alpha,Ap)
        #rel_err = norm(r)/normr0
        rel_err = Sparso.norm(r)/normr0
        blas1_time += time()

        if rel_err < tol 
            break
        end

        blas1_time -= time()
        #z = copy(r)
        Sparso.copy!(z, r)
        blas1_time += time()

        trsv_time -= time()
        #Base.SparseMatrix.fwdTriSolve!(L, z)
        Sparso.fwdTriSolve!(L, z, __fknob_8201)
        trsv_time += time()

        Sparso.reordering(
            __fknob_8201, 
            reordering_status, 
            A, Sparso.COL_PERM, Sparso.ROW_INV_PERM, __mknobA,
            U, Sparso.COL_PERM, Sparso.ROW_INV_PERM, __mknobU,
            :__delimitor__,
            r, Sparso.ROW_PERM,
            x, Sparso.COL_PERM,
            p, Sparso.COL_PERM
        )

        trsv_time -= time()
        #Base.SparseMatrix.bwdTriSolve!(U, z)
        Sparso.bwdTriSolve!(U, z, __fknob_8221)
        trsv_time += time()

        blas1_time -= time()
        #rz = dot(r, z)
        rz = Sparso.dot(r,z)
        beta = rz/old_rz
        #p = z + beta * p
        Sparso.WAXPBY!(p,1,z,beta,p)
        blas1_time += time()

        k += 1
    end

    Sparso.reverse_reordering(
        reordering_status, 
        :__delimitor__, 
        x, Sparso.COL_PERM
    )

    (Sparso.delete_function_knob)(__fknob_8221)
    (Sparso.delete_function_knob)(__fknob_8201)
    (Sparso.delete_matrix_knob)(__mknobL)

    if do_print
      total_time = time() - total_time
      println("total = $(total_time)s trsv_time = $(trsv_time)s ($((12.*(nnz(L) + nnz(U)) + 2.*8*(size(L, 1) + size(L, 2)))*k/trsv_time/1e9) gbps) spmv_time = $(spmv_time)s ($((12.*nnz(A) + 8.*(size(A, 1) + size(A, 2)))*(k + 1)/spmv_time/1e9) gbps) blas1_time = $blas1_time reorder_time = $reorder_time")
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

x, k, rel_err = pcg_symgs_ilu0_with_context_opt_without_reordering(x, A, b, tol, maxiter, false)
x = zeros(m)
println("Opt: ")
x, k, rel_err = pcg_symgs_ilu0_with_context_opt_without_reordering(x, A, b, tol, maxiter, true)
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

x, k, rel_err = pcg_symgs_ilu0_with_context_opt(x, A, b, tol, maxiter, false)
x = zeros(m)
println("Opt_with_reordering: ")
#Sparso.set_knob_log_level(1)
x, k, rel_err = pcg_symgs_ilu0_with_context_opt(x, A, b, tol, maxiter, true)
#Sparso.set_knob_log_level(0)
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

A2 = copy(A) # workaround that we change A in-place
@acc x, k, rel_err = pcg_symgs_ilu0(x, A2, b, tol, maxiter)
println("\nAccelerated: ")
x = zeros(m)
#Sparso.set_knob_log_level(1)
@acc x, k, rel_err = pcg_symgs_ilu0(x, A, b, tol, maxiter)
println("\tAccelerated k=", k)
println("\tAccelerated rel_err=", rel_err)
