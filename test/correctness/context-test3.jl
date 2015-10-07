include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS)

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
SparseAccelerator.set_reuse_inspection(false)
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
SparseAccelerator.set_knob_log_level(1)
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
@printf "\nAccelerated acc_total_time = %f\n" ref_total_time
@printf "Accelerated spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "Accelerated iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval
SparseAccelerator.set_reuse_inspection(true)

println("\n\nAccelerated: ")
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
SparseAccelerator.set_knob_log_level(1)
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)
@printf "\nAccelerated acc_total_time = %f\n" ref_total_time
@printf "Accelerated spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "Accelerated iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval

#println("\n\nWith manual context-sensitive optimization: ")
#ipm_ref_simplified_with_context_opt(A, b, p) 
#println("\tsum of x=", sum(x))
