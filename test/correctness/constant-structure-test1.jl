include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("./ipm-ref.jl")

A = matrix_market_read(string(ARGS[1], "-A.mtx"))'
b = vec(matrix_market_read(string(ARGS[1], "-b.mtx")))
p = vec(matrix_market_read(string(ARGS[1], "-p.mtx")))

println("\n\nAccelerated: ")
@acc x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_ref(A, b, p)

