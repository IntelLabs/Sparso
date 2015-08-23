function speye_int32(m::Integer)
  rowval = [Int32(x) for x in [1:m;]]
  colptr = [Int32(x) for x in [rowval; m + 1]]
  nzval = ones(Float64, m)
  return SparseMatrixCSC(m, m, colptr, rowval, nzval)
end


# a reference implementation of interior-point method
# that uses Julia's default SpGEMM and Cholesky
function ipm_ref1(A, b, p, print = true) # A: constraint coefficients, b: constraint rhs, p: objective
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

function ipm_ref(A, b, p) # A: constraint coefficients, b: constraint rhs, p: objective
  (m,n) = size(A)

  # set initial point, based on largest element in (A,b,p)
  bigM = maximum(A)
  bigM = maximum([norm(b, Inf) norm(p, Inf) bigM])
  x = 100*bigM*ones(n)
  s = x
  y = zeros(m)
  bc = 1 + maximum([norm(b) norm(p)])
  relResidual = NaN
  D = speye_int32(n)
  for iter=1:200
    # compute residuals
    Rd = A'*y + s - p
    Rp = A*x - b
    Rc = x.*s
    mu = mean(Rc)
    relResidual = norm([Rd; Rp; Rc])/bc
    if (relResidual <= 1e-7 && mu <= 1e-7) break; end

    Rc = Rc - min(0.1, 100*mu)*mu

    # set up the scaling matrix, and form the coefficient matrix for
    # the linear system
    d = min(5.e+15, x./s)

    D.nzval = d
    B = A*D*A'

    # use the form of the Cholesky routine "cholinc" that's best
    # suited to interior-point methods
    R = cholfact(SparseMatrixCSC{Float64, Int64}(B))

    # set up the right-hand side
    t1 = x.*Rd - Rc;
    t2 = -(Rp + A*(t1./s));

    # solve it and recover the other step components
    dy = R\t2
    temp = A'*dy
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
  end

  f = p'*x

  x
end