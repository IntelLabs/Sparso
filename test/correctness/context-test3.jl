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

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS)#, SA_DISABLE_COLLECTIVE_STRUCTURE_PREDICTION)

include("./ipm-ref.jl")
include("utils.jl")

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
@printf "\nOriginal ref_total_time = %f\n" ref_total_time
@printf "Original spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "Original iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval

println("\n\nWithout context-sensitive optimization: ")
Sparso.set_reuse_inspection(false)
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
Sparso.set_knob_log_level(1)
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
@printf "\nAccelerated acc_total_time = %f\n" ref_total_time
@printf "Accelerated spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "Accelerated iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval
Sparso.set_reuse_inspection(true)

println("\n\nAccelerated: ")
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
Sparso.set_knob_log_level(1)
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
@printf "\nAccelerated acc_total_time = %f\n" ref_total_time
@printf "Accelerated spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "Accelerated iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval

#println("\n\nWith manual context-sensitive optimization: ")
#ipm_ref_simplified_with_context_opt(A, b, p) 
#println("\tsum of x=", sum(x))
