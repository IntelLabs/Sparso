include("../../src/SparseAccelerator.jl")
using SparseAccelerator

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
  __mknob__AT___8599 = SparseAccelerator.new_matrix_knob(__AT__,true,true,false,false,false,false)
  __mknobA_8601      = SparseAccelerator.new_matrix_knob(A,     true,true,false,false,false,false)

  # Create mknobs for constant-structured matrices in the source code
  # Note R and B are also single-defs: key to enable context info propagation.
  # D is also single-def, but that is not important for this application.
  __mknobR_8596 = SparseAccelerator.new_matrix_knob(false,true,false,false,false,true)
  __mknobB_8598 = SparseAccelerator.new_matrix_knob(false,true,false,false,false,true)
  __mknobD_8600 = SparseAccelerator.new_matrix_knob(false,true,false,false,false,true)

  # Create mknobs for a function call at each call site, representing the
  # call's output, if the output is constant in structure or value. We do not
  # neeed this for cholmod_factor_inverse_divide, since its output is a vector,
  # while we track only matrices: it is OK, because our principle is: if there
  # is context info, make sure it is always correct and we can use it; if there
  # is not, we simply use nothing. Any way, a function call's result must remain
  # the same, no matter it has a fknob with it or not.
  # Any function call's result is a single-def (Defined only the function call).
  __mknob_ADB            = SparseAccelerator.new_matrix_knob(false,true,false,false,false,true)
  __mknob_cholfact_int32 = SparseAccelerator.new_matrix_knob(false,true,false,false,false,true)

  # Create fknobs for a function call at each call site, and add related mknobs
  # For ADB
  __fknob_8602 = SparseAccelerator.new_function_knob()
  (SparseAccelerator.add_mknob_to_fknob)(__mknob_ADB,__fknob_8602)
  (SparseAccelerator.add_mknob_to_fknob)(__mknob__AT___8599,__fknob_8602)
  (SparseAccelerator.add_mknob_to_fknob)(__mknobD_8600,__fknob_8602)
  (SparseAccelerator.add_mknob_to_fknob)(__mknobA_8601,__fknob_8602)

  (SparseAccelerator.set_derivative)(__mknobA_8601, SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE, __mknob__AT___8599)

  # For cholfact_int32
  __fknob_8622 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknobB_8598,__fknob_8622)

  # For cholmod_factor_inverse_divide
  __fknob_8623 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknobR_8596,__fknob_8623)

  fknob_spmv1 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknob__AT___8599, fknob_spmv1)

  fknob_spmv2 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknobA_8601, fknob_spmv2)

  fknob_spmv3 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknobA_8601, fknob_spmv3)

  fknob_spmv4 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknob__AT___8599, fknob_spmv4)

  for iter=1:200

    # compute residuals
    spmv_time -= time()
    #Rd = __AT__*y + s - p
    Rd = SparseAccelerator.SpMV(1, __AT__, y, 1, s - p, fknob_spmv1)
    #Rp = A*x - b
    Rp = SparseAccelerator.SpMV(1, A, x, -1, b, fknob_spmv2)
    spmv_time += time()

    blas1_time -= time()
    #Rc = x.*s
    Rc = SparseAccelerator.element_wise_multiply(x, s)
    #mu = mean(Rc)
    mu = SparseAccelerator.sum(Rc)/length(Rc)
    relResidual = norm([Rd; Rp; Rc])/bc
    blas1_time += time()

    if (relResidual <= 1e-7 && mu <= 1e-7) break; end

    blas1_time -= time()
    #Rc = Rc - min(0.1, 100*mu)*mu
    SparseAccelerator.WAXPB!(Rc, 1, Rc, -min(0.1, 100*mu)*mu)

    # set up the scaling matrix, and form the coefficient matrix for
    # the linear system
    #d = min(5.e+15, x./s)
    d = SparseAccelerator.element_wise_divide(x, s)
    SparseAccelerator.min!(d, d, 5.e+15)
    blas1_time += time()

    spgemm_time -= time()
    D.nzval = d
    #B = A*D*A'
    B = SparseAccelerator.ADB(A,D,__AT__,__fknob_8602)
    SparseAccelerator.propagate_matrix_info(__mknobB_8598, __mknob_ADB)
    spgemm_time += time()

    # use the form of the Cholesky routine "cholinc" that's best
    # suited to interior-point methods
    fact_time -= time()
    #R = cholfact_int32(B)
    R = SparseAccelerator.cholfact_int32(B,__fknob_8622)
    fact_time += time()

    # set up the right-hand side
    blas1_time -= time()
    #t1 = x.*Rd - Rc;
    t1 = SparseAccelerator.element_wise_multiply(x, Rd)
    SparseAccelerator.WAXPBY!(t1, 1, t1, -1, Rc)
    blas1_time += time()

    spmv_time -= time()
    #t2 = -(Rp + A*(t1./s));
    t2 = SparseAccelerator.SpMV(-1, A, t1./s, -1, Rp, fknob_spmv3)
    spmv_time += time()

    # solve it and recover the other step components
    trslv_time -= time()
    #dy = R\t2
    dy = SparseAccelerator.cholfact_inverse_divide(R,t2,__fknob_8623)
    trslv_time += time()

    spmv_time -= time()
    #temp = A'*dy
    temp = SparseAccelerator.SpMV(__AT__, dy, fknob_spmv4)
    spmv_time += time()

    blas1_time -= time()
    #dx = (x.*temp + t1)./s
    SparseAccelerator.element_wise_multiply!(temp, x, temp)
    SparseAccelerator.WAXPBY!(temp, 1, temp, 1, t1)
    dx = SparseAccelerator.element_wise_divide(temp, s)
    #ds = -(s.*dx + Rc)./x
    SparseAccelerator.element_wise_multiply!(temp, s, dx)
    SparseAccelerator.WAXPBY!(temp, -1, temp, -1, Rc)
    ds = SparseAccelerator.element_wise_divide(temp, x)

    tau = max(.9995, 1 - mu)
    #ap = -1/minimum([dx./x; -1])
    SparseAccelerator.element_wise_divide!(temp, dx, x)
    ap = -1/min(SparseAccelerator.minimum(temp), -1)
    #ad = -1/minimum([ds./s; -1])
    SparseAccelerator.element_wise_divide!(temp, ds, s)
    ad = -1/min(SparseAccelerator.minimum(temp), -1)
    ap = tau*ap
    ad = tau*ad
    #x = x + ap*dx
    #s = s + ad*ds
    #y = y + ad*dy
    SparseAccelerator.WAXPBY!(x, 1, x, ap, dx)
    SparseAccelerator.WAXPBY!(s, 1, s, ad, ds)
    SparseAccelerator.WAXPBY!(y, 1, y, ad, dy)
    blas1_time += time()
  end

  (SparseAccelerator.delete_function_knob)(__fknob_8602)
  (SparseAccelerator.delete_function_knob)(__fknob_8622)
  (SparseAccelerator.delete_function_knob)(__fknob_8623)
  (SparseAccelerator.delete_matrix_knob)(__mknobR_8596)
  (SparseAccelerator.delete_matrix_knob)(__mknobB_8598)
  (SparseAccelerator.delete_matrix_knob)(__mknob__AT___8599)
  (SparseAccelerator.delete_matrix_knob)(__mknobA_8601)
  (SparseAccelerator.delete_matrix_knob)(__mknobD_8600) # line 56:
  (SparseAccelerator.delete_matrix_knob)(__mknob_ADB)
  (SparseAccelerator.delete_matrix_knob)(__mknob_cholfact_int32) # line 56:

  blas1_time -= time()
  f = p'*x
  blas1_time += time()

  x, time() - t0, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time, iter, relResidual, f[1][1]
end

if length(ARGS) == 0
    println("Need an input matrix. Try add ipm/mps/osa-14")
    assert(false)
    # This does not seem to work.
    # min 2*x1 + x2 subject to x1 + x2 = 1, x1 >= 0, x2 >= 0
    # expected solution: x1 = 0, x2 = 1, obj = 1
    A = sparse([1 1])
    b = [ 1 ]'
    p = [ 2 1 ]'
else
    A = matrix_market_read(string(ARGS[1], "-A.mtx"))'
    b = matrix_market_read(string(ARGS[1], "-b.mtx"))
    p = matrix_market_read(string(ARGS[1], "-p.mtx"))
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
println("Original sum of x=", sum(x))
@printf "\nref_total_time = %f\n" ref_total_time
@printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval

println("\n\nWith manual context-sensitive optimization: ")
x, ref_total_time, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time,
    iter, relResidual, objval = ipm_with_context_opt(A, b, p)
println("Manual_context sum of x=", sum(x))
@printf "\nopt_total_time = %f\n" ref_total_time
@printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
@printf "iter %2i, resid = %9.2e, objval = %e\n" iter relResidual objval

#println("\n\nWith manual context-sensitive optimization: ")
#ipm_ref_simplified_with_context_opt(A, b, p) 
#println("\tsum of x=", sum(x))



