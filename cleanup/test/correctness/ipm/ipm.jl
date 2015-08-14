using MatrixMarket
using MatrixMarket2

include("dss.jl")

function speye_int32(m::Integer)
  rowval = [Int32(x) for x in [1:m;]]
  colptr = [Int32(x) for x in [rowval; m + 1]]
  nzval = ones(Float64, m)
  return SparseMatrixCSC(m, m, colptr, rowval, nzval)
end

# A*D*B, where D is a diagonal matrix
# Since it's implemented in CSR, we pass tranposes of A and B
# in CSC
function adb_inspect(AT::SparseMatrixCSC, BT::SparseMatrixCSC)
  csrA = SparseAccelerator.CreateCSR(AT)
  csrB = SparseAccelerator.CreateCSR(BT)

  # A*D*A' inspection
  csrADAT = ccall((:CSR_ADBInspect, LIB_CSR_PATH), Ptr{Void},
                  (Ptr{Void}, Ptr{Void}),
                  csrA, csrB)

  m = size(A, 1)
  rowptr = pointer_to_array(
    ccall((:CSR_GetRowPtr, LIB_CSR_PATH), Ptr{Cint},
          (Ptr{Void},),
          csrADAT),
    (m + 1,))
  nnz = rowptr[m + 1] - 1
  colidx = pointer_to_array(
    ccall((:CSR_GetColIdx, LIB_CSR_PATH), Ptr{Cint},
          (Ptr{Void},),
          csrADAT),
    (nnz,))
  values = pointer_to_array(
    ccall((:CSR_GetValues, LIB_CSR_PATH), Ptr{Cdouble},
          (Ptr{Void},),
          csrADAT),
    (nnz,))

  SparseAccelerator.DestroyCSR(csrA)
  SparseAccelerator.DestroyCSR(csrB)

  ADB = SparseMatrixCSC{Cdouble, Cint}(
    m, m, rowptr, colidx, values)
end

# A*D*B, where D is a diagonal matrix
function adb_execute!(ADB::SparseMatrixCSC, AT::SparseMatrixCSC, BT::SparseMatrixCSC, d::Vector)
  csrA = SparseAccelerator.CreateCSR(AT)
  csrB = SparseAccelerator.CreateCSR(BT)
  csrADB = SparseAccelerator.CreateCSR(ADB)

  ccall((:CSR_ADB, LIB_CSR_PATH), Void,
        (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Cdouble}),
        csrADB, csrA, csrB, d)

  SparseAccelerator.DestroyCSR(csrA)
  SparseAccelerator.DestroyCSR(csrB)
  SparseAccelerator.DestroyCSR(csrADB)
end

# optimized implementation of interior-point method with inspector hoisted
function ipm(A, b, p, print = true) # A: constraint coefficients, b: constraint rhs, p: objective
  t0 = time()
  (m,n) = size(A)

# set initial point, based on largest element in (A,b,p)
  bigM = maximum(A)
  bigM = maximum([norm(b, Inf) norm(p, Inf) bigM])
  x = 100*bigM*ones(n)
  s = x
  y = zeros(m)

  bc = 1 + maximum([norm(b) norm(p)])

  relResidual = NaN
  iter=1

  blas1_time = 0.
  spgemm_time = 0.
  fact_time = 0.
  trslv_time = 0.
  spmv_time = 0.

  spmv_inspector_time = 0.
  spgemm_inspector_time = 0.
  fact_inspector_time = 0.

  D = speye_int32(n)

  spmv_inspector_time -= time()
  AT = A'
  spmv_inspector_time += time()

  spgemm_inspector_time -= time()
  ADAT = adb_inspect(AT, A)
  spgemm_inspector_time += time()

  fact_inspector_time -= time()

  dss_handle = dss_analyze(ADAT)
  dy = Array(Cdouble, m)

  fact_inspector_time += time()

  for iter=1:200

    # compute residuals
    spmv_time -= time()
    Rd = AT*y + s - p
    Rp = A*x - b
    spmv_time += time()

    blas1_time -= time()
    Rc = x.*s
    mu = mean(Rc)
    relResidual = norm([Rd; Rp; Rc])/bc
    blas1_time += time()

    if print
      #@printf "iter %2i: log10(mu) = %5.2f, resid = %9.2e, obj = %9.2e\n" iter log10(mu) relResidual (p'*x)[1][1]
    end
    if (relResidual <= 1e-7 && mu <= 1e-7) break; end

    blas1_time -= time()
    Rc = Rc - min(0.1, 100*mu)*mu

    # set up the scaling matrix, and form the coefficient matrix for
    # the linear system
    d = min(5.e+15, x./s)
    blas1_time += time()

    spgemm_time -= time()
    adb_execute!(ADAT, AT, A, d)
    spgemm_time += time()

    # use the form of the Cholesky routine "cholinc" that's best
    # suited to interior-point methods
    fact_time -= time()
    dss_factor(dss_handle, ADAT)
    fact_time += time()

    # set up the right-hand side
    blas1_time -= time()
    t1 = x.*Rd - Rc;
    blas1_time += time()

    spmv_time -= time()
    t2 = -(Rp + A*(t1./s));
    spmv_time += time()

    # solve it and recover the other step components
    trslv_time -= time()
    opt = MKL_DSS_DEFAULTS
    dss_solve!(dss_handle, t2, dy)
    trslv_time += time()

    spmv_time -= time()
    temp = AT*dy
    spmv_time += time()

    blas1_time -= time()
    dx = (x.*temp + t1)./s
    ds = -(s.*dx + Rc)./x

    tau = max(.9995, 1 - mu)
    ap = -1/minimum([dx./x; -1])
    ad = -1/minimum([ds./s; -1])
    ap = tau*ap
    ad = tau*ad
    x = x + ap*dx
    s = s + ad*ds
    y = y + ad*dy
    blas1_time += time()
  end

  blas1_time -= time()
  f = p'*x
  blas1_time += time()

  if print
    @printf "\nwith_inspector_hoisting_total_time = %f\n" time() - t0
    @printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
    @printf "spgemm_inspector = %f fact_inspector_time = %f spmv_inspector = %f\n" spgemm_inspector_time fact_inspector_time spmv_inspector_time
    @printf "iter %2i, resid = %9.2e, objval = %e\n" iter relResidual f[1][1]
  end

  x
end

# unoptimized implementation of interior-point method w/o inspector hoisting
function ipm_unopt(A, b, p, print = true) # A: constraint coefficients, b: constraint rhs, p: objective
  t0 = time()
  (m,n) = size(A)

# set initial point, based on largest element in (A,b,p)
  bigM = maximum(A)
  bigM = maximum([norm(b, Inf) norm(p, Inf) bigM])
  x = 100*bigM*ones(n)
  s = x
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
  dy = Array(Cdouble, m)

  for iter=1:200

    # compute residuals
    spmv_time -= time()
    Rd = A'*y + s - p
    Rp = A*x - b
    spmv_time += time()

    blas1_time -= time()
    Rc = x.*s
    mu = mean(Rc)
    relResidual = norm([Rd; Rp; Rc])/bc
    blas1_time += time()

    if print
      #@printf "iter %2i: log10(mu) = %5.2f, resid = %9.2e, obj = %9.2e\n" iter log10(mu) relResidual (p'*x)[1][1]
    end
    if (relResidual <= 1e-7 && mu <= 1e-7) break; end

    blas1_time -= time()
    Rc = Rc - min(0.1, 100*mu)*mu

    # set up the scaling matrix, and form the coefficient matrix for
    # the linear system
    d = min(5.e+15, x./s)
    blas1_time += time()

    spgemm_time -= time()
    ADAT = adb_inspect(A', A)
    adb_execute!(ADAT, A', A, d)
    spgemm_time += time()

    # use the form of the Cholesky routine "cholinc" that's best
    # suited to interior-point methods
    fact_time -= time()
    dss_handle = dss_analyze(ADAT)
    dss_factor(dss_handle, ADAT)
    fact_time += time()

    # set up the right-hand side
    blas1_time -= time()
    t1 = x.*Rd - Rc;
    blas1_time += time()

    spmv_time -= time()
    t2 = -(Rp + A*(t1./s));
    spmv_time += time()

    # solve it and recover the other step components
    trslv_time -= time()
    opt = MKL_DSS_DEFAULTS
    dss_solve!(dss_handle, t2, dy)
    trslv_time += time()

    spmv_time -= time()
    temp = A'*dy
    spmv_time += time()

    blas1_time -= time()
    dx = (x.*temp + t1)./s
    ds = -(s.*dx + Rc)./x

    tau = max(.9995, 1 - mu)
    ap = -1/minimum([dx./x; -1])
    ad = -1/minimum([ds./s; -1])
    ap = tau*ap
    ad = tau*ad
    x = x + ap*dx
    s = s + ad*ds
    y = y + ad*dy
    blas1_time += time()
  end

  blas1_time -= time()
  f = p'*x
  blas1_time += time()

  if print
    @printf "\nwo_inspector_hoisting_total_time = %f\n" time() - t0
    @printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
    @printf "iter %2i, resid = %9.2e, objval = %e\n" iter relResidual f[1][1]
  end

  x
end

# a reference implementation of interior-point method
# that uses Julia's default SpGEMM and Cholesky
function ipm_ref(A, b, p, print = true) # A: constraint coefficients, b: constraint rhs, p: objective
  t0 = time()
  (m,n) = size(A)

# set initial point, based on largest element in (A,b,p)
  bigM = maximum(A)
  bigM = maximum([norm(b, Inf) norm(p, Inf) bigM])
  x = 100*bigM*ones(n)
  s = x
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

  for iter=1:200

    # compute residuals
    spmv_time -= time()
    Rd = A'*y + s - p
    Rp = A*x - b
    spmv_time += time()

    blas1_time -= time()
    Rc = x.*s
    mu = mean(Rc)
    relResidual = norm([Rd; Rp; Rc])/bc
    blas1_time += time()

    if print
      #@printf "iter %2i: log10(mu) = %5.2f, resid = %9.2e, obj = %9.2e\n" iter log10(mu) relResidual (p'*x)[1][1]
    end
    if (relResidual <= 1e-7 && mu <= 1e-7) break; end

    blas1_time -= time()
    Rc = Rc - min(0.1, 100*mu)*mu

    # set up the scaling matrix, and form the coefficient matrix for
    # the linear system
    d = min(5.e+15, x./s)
    blas1_time += time()

    spgemm_time -= time()
    D.nzval = d
    B = A*D*A'
    spgemm_time += time()

    # use the form of the Cholesky routine "cholinc" that's best
    # suited to interior-point methods
    fact_time -= time()
    R = cholfact(SparseMatrixCSC{Float64, Int64}(B))
    fact_time += time()

    # set up the right-hand side
    blas1_time -= time()
    t1 = x.*Rd - Rc;
    blas1_time += time()

    spmv_time -= time()
    t2 = -(Rp + A*(t1./s));
    spmv_time += time()

    # solve it and recover the other step components
    trslv_time -= time()
    dy = R\t2
    trslv_time += time()

    spmv_time -= time()
    temp = A'*dy
    spmv_time += time()

    blas1_time -= time()
    dx = (x.*temp + t1)./s
    ds = -(s.*dx + Rc)./x

    tau = max(.9995, 1 - mu)
    ap = -1/minimum([dx./x; -1])
    ad = -1/minimum([ds./s; -1])
    ap = tau*ap
    ad = tau*ad
    x = x + ap*dx
    s = s + ad*ds
    y = y + ad*dy
    blas1_time += time()
  end

  blas1_time -= time()
  f = p'*x
  blas1_time += time()

  if print
    @printf "\nref_total_time = %f\n" time() - t0
    @printf "spgemm = %f fact = %f blas1 = %f trslv = %f spmv = %f\n" spgemm_time fact_time blas1_time trslv_time spmv_time
    @printf "iter %2i, resid = %9.2e, objval = %e\n" iter relResidual f[1][1]
  end

  x
end

# min 2*x1 + x2 subject to x1 + x2 = 1, x1 >= 0, x2 >= 0
# expected solution: x1 = 0, x2 = 1, obj = 1
#A = sparse([1 1])
#b = [ 1 ]'
#p = [ 2 1 ]'
A = MatrixMarket2.mmread(string(ARGS[1], "-A.mtx"))'
b = vec(MatrixMarket.mmread(string(ARGS[1], "-b.mtx")))
p = vec(MatrixMarket.mmread(string(ARGS[1], "-p.mtx")))

m = size(A, 1)
n = size(A, 2)
println("Problem size = [$m $n]")

ipm_ref(A, b, p, false) # ignore timing of the first run
ipm_ref(A, b, p)

ipm_unopt(A, b, p, false)
ipm_unopt(A, b, p)

ipm(A, b, p, false)
ipm(A, b, p)
