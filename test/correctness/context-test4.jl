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

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("./ipm-ref.jl")
include("utils.jl")

# This code is what we expect to generate by context-sensitive optimizations.
# It is used for debugging purpose only
function ipm_with_context_opt(A, b, p) # A: constraint coefficients, b: constraint rhs, p: objective
  t0 = time()
  (m,n) = size(A)

# set initial point, based on largest element in (A,b,p)
  bigM = maximum(A)
  bigM = maximum([norm(b, Inf) norm(p, Inf) bigM])
  x = 100*bigM*ones(n)
  s = copy(x)
  y = zeros(m)

  bc = 1 + maximum([norm(b) norm(p)])

  relResidual = NaN
  iter=1

  blas1_time = 0.
  spgemm_time = 0.
  fact_time = 0.
  trslv_time = 0.
  spmv_time = 0.

  D = speye_int32(n)
  
  # There are two kinds of matrices we would like to trace: one is constant value,
  # the other constant structure.
  # For a constant-valued matrix, since it is not defined but referenced in the
  # loop, it must already have a reprsentation in a fixed memory. Thus we can 
  # safely pass its name into new_matrix_knob(), which will remember its 
  # representation. 
  # For a constant-structured matrix, however, it is defined in the loop, and in
  # each definition, it might be assigned a representation in a different memory,
  # even though every time the representation assigned has the same structure.
  # Therefore, we cannot pass its name into new_matrix_knob() at this moment, and
  # new_matrix_knob() should not remember its representation, which may not 
  # even exist at this moment. Instead, every time is it defined in the loop, its
  # representation should be recorded into its matrix knob.

  # Hoist AT out of loop. This is currently done in pattern replacement.
  # TODO: have a separate loop invariant hoisting phase.
  __AT__ = (Main.ctranspose)(A::Base.SparseMatrix.SparseMatrixCSC{Float64,Int32})

  # Create mknobs for constant-valued matrices. They are of course constant-structured.
  __mknob__AT___8599 = Sparso.new_matrix_knob(:__AT__,true,true,false,false,false,false)
  __mknobA_8601      = Sparso.new_matrix_knob(:A,     true,true,false,false,false,false)

  # Create mknobs for constant-structured matrices in the source code
  # Note R and B are also single-defs: key to enable context info propagation.
  # D is also single-def, but that is not important for this application.
  __mknobR_8596 = Sparso.new_matrix_knob(:R, false,true,false,false,false,true)
  __mknobB_8598 = Sparso.new_matrix_knob(:B, false,true,false,false,false,true)
  __mknobD_8600 = Sparso.new_matrix_knob(:D, false,true,false,false,false,true)

  # Create mknobs for a function call at each call site, representing the
  # call's output, if the output is constant in structure or value. We do not
  # neeed this for cholmod_factor_inverse_divide, since its output is a vector,
  # while we track only matrices: it is OK, because our principle is: if there
  # is context info, make sure it is always correct and we can use it; if there
  # is not, we simply use nothing. Any way, a function call's result must remain
  # the same, no matter it has a fknob with it or not.
  # Any function call's result is a single-def (Defined only the function call).
  __mknob_ADB            = Sparso.new_matrix_knob(:ADB, false,true,false,false,false,true)
  __mknob_cholfact_int32 = Sparso.new_matrix_knob(:cholfact_int32, false,true,false,false,false,true)

  # Create fknobs for a function call at each call site, and add related mknobs
  # For ADB
  __fknob_8602 = Sparso.new_function_knob()
  (Sparso.add_mknob_to_fknob)(__mknob_ADB,__fknob_8602)
  (Sparso.add_mknob_to_fknob)(__mknob__AT___8599,__fknob_8602)
  (Sparso.add_mknob_to_fknob)(__mknobD_8600,__fknob_8602)
  (Sparso.add_mknob_to_fknob)(__mknobA_8601,__fknob_8602)

  (Sparso.set_derivative)(__mknobA_8601, Sparso.DERIVATIVE_TYPE_TRANSPOSE, __mknob__AT___8599)

  # For cholfact_int32
  __fknob_8622 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobB_8598,__fknob_8622)

  # For cholmod_factor_inverse_divide
  __fknob_8623 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobR_8596,__fknob_8623)

  fknob_spmv1 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknob__AT___8599, fknob_spmv1)

  fknob_spmv2 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobA_8601, fknob_spmv2)

  fknob_spmv3 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobA_8601, fknob_spmv3)

  fknob_spmv4 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknob__AT___8599, fknob_spmv4)

  for iter=1:200

    # compute residuals
    spmv_time -= time()
    #Rd = __AT__*y + s - p
    Rd = Sparso.SpMV(1, __AT__, y, 1, s - p, fknob_spmv1)
    #Rp = A*x - b
    Rp = Sparso.SpMV(1, A, x, -1, b, fknob_spmv2)
    spmv_time += time()

    blas1_time -= time()
    #Rc = x.*s
    Rc = Sparso.element_wise_multiply(x, s)
    #mu = mean(Rc)
    mu = Sparso.sum(Rc)/length(Rc)
    relResidual = norm([Rd; Rp; Rc])/bc
    blas1_time += time()

    if (relResidual <= 1e-7 && mu <= 1e-7) break; end

    blas1_time -= time()
    #Rc = Rc - min(0.1, 100*mu)*mu
    Sparso.WAXPB!(Rc, 1, Rc, -min(0.1, 100*mu)*mu)

    # set up the scaling matrix, and form the coefficient matrix for
    # the linear system
    #d = min(5.e+15, x./s)
    d = Sparso.element_wise_divide(x, s)
    Sparso.min!(d, d, 5.e+15)
    blas1_time += time()

    spgemm_time -= time()
    D.nzval = d
    #B = A*D*A'
    B = Sparso.ADB(A,D,__AT__,__fknob_8602)
    Sparso.propagate_matrix_info(__mknobB_8598, __mknob_ADB)
    spgemm_time += time()

    # use the form of the Cholesky routine "cholinc" that's best
    # suited to interior-point methods
    fact_time -= time()
    #R = cholfact_int32(B)
    R = Sparso.cholfact_int32(B,__fknob_8622)
    fact_time += time()

    # set up the right-hand side
    blas1_time -= time()
    #t1 = x.*Rd - Rc;
    t1 = Sparso.element_wise_multiply(x, Rd)
    Sparso.WAXPBY!(t1, 1, t1, -1, Rc)
    blas1_time += time()

    spmv_time -= time()
    #t2 = -(Rp + A*(t1./s));
    t2 = Sparso.SpMV(-1, A, t1./s, -1, Rp, fknob_spmv3)
    spmv_time += time()

    # solve it and recover the other step components
    trslv_time -= time()
    #dy = R\t2
    dy = Sparso.cholfact_inverse_divide(R,t2,__fknob_8623)
    trslv_time += time()

    spmv_time -= time()
    #temp = A'*dy
    temp = Sparso.SpMV(__AT__, dy, fknob_spmv4)
    spmv_time += time()

    blas1_time -= time()
    #dx = (x.*temp + t1)./s
    Sparso.element_wise_multiply!(temp, x, temp)
    Sparso.WAXPBY!(temp, 1, temp, 1, t1)
    dx = Sparso.element_wise_divide(temp, s)
    #ds = -(s.*dx + Rc)./x
    Sparso.element_wise_multiply!(temp, s, dx)
    Sparso.WAXPBY!(temp, -1, temp, -1, Rc)
    ds = Sparso.element_wise_divide(temp, x)

    tau = max(.9995, 1 - mu)
    #ap = -1/minimum([dx./x; -1])
    Sparso.element_wise_divide!(temp, dx, x)
    ap = -1/min(Sparso.minimum(temp), -1)
    #ad = -1/minimum([ds./s; -1])
    Sparso.element_wise_divide!(temp, ds, s)
    ad = -1/min(Sparso.minimum(temp), -1)
    ap = tau*ap
    ad = tau*ad
    #x = x + ap*dx
    #s = s + ad*ds
    #y = y + ad*dy
    Sparso.WAXPBY!(x, 1, x, ap, dx)
    Sparso.WAXPBY!(s, 1, s, ad, ds)
    Sparso.WAXPBY!(y, 1, y, ad, dy)
    blas1_time += time()
  end

  (Sparso.delete_function_knob)(__fknob_8602)
  (Sparso.delete_function_knob)(__fknob_8622)
  (Sparso.delete_function_knob)(__fknob_8623)
  (Sparso.delete_matrix_knob)(__mknobR_8596)
  (Sparso.delete_matrix_knob)(__mknobB_8598)
  (Sparso.delete_matrix_knob)(__mknob__AT___8599)
  (Sparso.delete_matrix_knob)(__mknobA_8601)
  (Sparso.delete_matrix_knob)(__mknobD_8600) # line 56:
  (Sparso.delete_matrix_knob)(__mknob_ADB)
  (Sparso.delete_matrix_knob)(__mknob_cholfact_int32) # line 56:

  blas1_time -= time()
  f = p'*x
  blas1_time += time()

  x, time() - t0, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time, iter, relResidual, f[1][1]
end

if length(ARGS) == 0
  A, b, p = load_ipm_input("ipm/mps/osa-14")
else
  A, b, p = load_ipm_input(ARGS[1])
end

m = size(A, 1)
n = size(A, 2)
println("Problem size = [$m $n]")
#println("\tsum of A=", sum(A))
#println("\tsum of b=", sum(b))
#println("\tsum of p=", sum(p))

println("Original: ")
x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
@printf "Original iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval
@printf "\nref_total_time = %f\n" ref_total_time
@printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time

println("\n\nWith manual context-sensitive optimization: ")
x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_with_context_opt(A, b, p)
Sparso.set_knob_log_level(1)
x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_with_context_opt(A, b, p)
@printf "\nopt_total_time = %f\n" ref_total_time
@printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "Manual_context iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval

#println("\n\nWith manual context-sensitive optimization: ")
#ipm_ref_simplified_with_context_opt(A, b, p) 
#println("\tsum of x=", sum(x))
