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

# julia knob-test2.jl ../matrices/bcsstk14.mtx
# Turning on any of the commented statement to repro the bug.
 
include("../../src/Sparso.jl")
include("../../src/simple-show.jl")
include("./utils.jl")
using Sparso

function pcg_symgs(x, A, b, tolerance, maxiter)
    __mknobL4__ = Sparso.new_matrix_knob(:L,true,true,false,false,false,true)
    __mknobA1__ = Sparso.new_matrix_knob(:A,true,true,true,true,false,true)
    __mknobU6__ = Sparso.new_matrix_knob(:U,true,true,false,false,false,true)
    __fknobilu__ = Sparso.new_function_knob()
    __fknob3__ = Sparso.new_function_knob()
    __fknob2__ = Sparso.new_function_knob()
    __fknob5__ = Sparso.new_function_knob()
    __fknob7__ = Sparso.new_function_knob()
    __fknob0__ = Sparso.new_function_knob()
    __fknob8__ = Sparso.new_function_knob()
    Sparso.set_derivative(__mknobL4__,Sparso.DERIVATIVE_TYPE_LOWER_TRIANGULAR,__mknobA1__)
    Sparso.set_derivative(__mknobU6__,Sparso.DERIVATIVE_TYPE_SYMMETRIC,__mknobA1__)
    Sparso.set_derivative(__mknobL4__,Sparso.DERIVATIVE_TYPE_SYMMETRIC,__mknobA1__)
    Sparso.set_derivative(__mknobU6__,Sparso.DERIVATIVE_TYPE_UPPER_TRIANGULAR,__mknobA1__)
    Sparso.add_mknob_to_fknob(__mknobL4__,__fknob3__)
    Sparso.add_mknob_to_fknob(__mknobA1__,__fknob2__)
    Sparso.add_mknob_to_fknob(__mknobU6__,__fknob5__)
    Sparso.add_mknob_to_fknob(__mknobL4__,__fknob7__)
    Sparso.add_mknob_to_fknob(__mknobA1__,__fknob0__)
    Sparso.add_mknob_to_fknob(__mknobU6__,__fknob8__)
    Sparso.add_mknob_to_fknob(__mknobA1__,__fknobilu__)
    total_time = time()
    
    
    L = tril(A)
    U = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)
    # ISSUE: Either of the following will make the code run endlessly
    # L, U = Sparso.ilu(A)
    # L, U = Sparso.ilu(A, __fknobilu__)
    
    r = b - Sparso.SpMV(1,A,x,__fknob2__)
    normr0 = Sparso.norm(r)
    z = copy(r)
    Sparso.fwdTriSolve!(z,L,r,__fknob3__)
    Sparso.bwdTriSolve!(z,U,z,__fknob5__)
    p = copy(z)
    rz = Sparso.dot(r,z)
    rel_err = 1
    k = 1
    while true
            old_rz = rz
            Ap = Sparso.SpMV(1,A,p,__fknob0__)
            alpha = old_rz / Sparso.dot(p,Ap)
            Sparso.WAXPBY!(x,1,x,alpha,p)
            Sparso.WAXPBY!(r,1,r,-alpha,Ap)
            rel_err = Sparso.norm(r) / normr0
            if rel_err < tolerance
                break
            end
            Sparso.fwdTriSolve!(z,L,r,__fknob7__)
            Sparso.bwdTriSolve!(z,U,z,__fknob8__)
            rz = Sparso.dot(r,z)
            Sparso.WAXPBY!(p,1,z,rz / old_rz,p)
            k = k + 1
    end
    total_time = time() - total_time
    println("total = ",total_time,"s")
    Sparso.delete_matrix_knob(__mknobL4__)
    Sparso.delete_matrix_knob(__mknobA1__)
    Sparso.delete_matrix_knob(__mknobU6__)
    Sparso.delete_function_knob(__fknob3__)
    Sparso.delete_function_knob(__fknob2__)
    Sparso.delete_function_knob(__fknob5__)
    Sparso.delete_function_knob(__fknob7__)
    Sparso.delete_function_knob(__fknob0__)
    Sparso.delete_function_knob(__fknob8__)
    return x,k,rel_err
end

A         = matrix_market_read(ARGS[1], true, true)
m         = size(A, 1)
b         = ones(Float64, m)
originalA = copy(A)
tolerance = 1e-7
maxiter   = 20000

x = zeros(Float64, m)
x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
println("\tOptimized sum of x=", sum(x))
println("\tOptimized k=", k)
println("\tOptimized rel_err=", rel_err)
