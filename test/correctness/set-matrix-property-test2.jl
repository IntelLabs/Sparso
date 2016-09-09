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
include("./utils.jl")
using Sparso

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

# The original pcg_symgs
#include("./pcg-symgs.jl")
function pcg_symgs(x, A, b, tol, maxiter)
    total_time = time()
  
    trsv_time = 0.
    spmv_time = 0.
    dot_time = 0.
    blas1_time = 0.

    L = tril(A)
    U  = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A))*triu(A))

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
    while k <= maxiter

        set_matrix_property(Dict(
            :A => SA_CONST_VALUED|SA_MAXIMAL_STRUCTURED|SA_SYMM_VALUED|SA_SYMM_STRUCTURED, 
            :L => SA_CONST_VALUED|SA_MAXIMAL_STRUCTURED,
            :U => SA_CONST_VALUED|SA_MAXIMAL_STRUCTURED
            )
        )

        old_rz = rz

        spmv_time -= time()
        #Ap = A*p
        Ap = Sparso.SpMV(A,p)
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
        Base.SparseMatrix.fwdTriSolve!(L, z)
        Base.SparseMatrix.bwdTriSolve!(U, z)
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
tol     = 1e-7
maxiter = 20000

println("Original: ")
x, k, rel_err = @acc pcg_symgs(x, A, b, tol, maxiter)
#println("\tsum of x=", sum(x))
println("\tk=", k)
println("\trel_err=", rel_err)


