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

function pcg_jacobi(x, A, b, tol, maxiter)
    inv_d = 1./diag(A)
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = inv_d .* r
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    while k <= maxiter
        old_rz = rz
        Ap = A*p #Ap = Sparso.SpMV(A, p) # manual # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap) #alpha = old_rz / Sparso.dot(p, Ap) # manual
        x += alpha * p #Sparso.WAXPBY!(x, alpha, p, 1, x) # manual
        r -= alpha * Ap #Sparso.WAXPBY!(r, -alpha, Ap, 1, r) # manual
        rel_err = sqrt(dot(r, r))/normr0 #rel_err = sqrt(Sparso.dot(r, r))/normr0 # manual
        if rel_err < tol 
            break
        end
        z = inv_d .* r #z = Sparso.element_wise_multiply(inv_d, r)
        rz = dot(r, z) #rz = Sparso.dot(r, z) # manual
        beta = rz/old_rz
        p = z + beta * p #Sparso.WAXPBY!(p, 1, z, beta, p) # manual
        k += 1
    end
    return x, k, rel_err
end
