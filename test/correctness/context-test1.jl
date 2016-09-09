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

#include("./pcg-symgs.jl")
include("utils.jl")

# TEMPORARY: this should be removed in future once structure analysis can
# recognize that L and U are part of A
# Copied from pcg-symgs.jl to avoid affecting other tests if it were changed.
function pcg_symgs(x, A, b, tol, maxiter)
    # These are the only 3 properties we need to discover automatically.
    # We also require SA_CONST_STRUCTURED for A, but that has already been discovered.
#    set_matrix_property(:L, SA_LOWER_OF, :A)
#    set_matrix_property(:U, SA_UPPER_OF, :A)
    set_matrix_property(Dict(
        :A => SA_SYMM_STRUCTURED | SA_SYMM_VALUED,
        :U => SA_CONST_VALUED
        )
    )

    total_time = time()

    trsv_time = 0.
    spmv_time = 0.
    blas1_time = 0.

    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A))*triu(A))

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
        Ap = A*p # Ap = Sparso.SpMV(A, p) # This takes most time. Compiler can reorder A to make faster
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
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-7
maxiter = 20000

x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
x = zeros(Float64, m)
println("Original: ")
x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

A2 = copy(A) # workaround that we change A in-place
@acc x, k, rel_err = pcg_symgs(x, A2, b, tol, maxiter)
println("\nAccelerated: ")
x = zeros(Float64, m)
@acc x, k, rel_err = pcg_symgs(x, A, b, tol, maxiter)
println("\tAccelerated k=", k)
println("\tAccelerated rel_err=", rel_err)
