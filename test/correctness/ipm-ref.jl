include("ipm/dss.jl")

function speye_int32(m::Integer)
  rowval = [Int32(x) for x in [1:m;]]
  colptr = [Int32(x) for x in [rowval; m + 1]]
  nzval = ones(Float64, m)
  return SparseMatrixCSC(m, m, colptr, rowval, nzval)
end

function canonicalize_ipm_input(A, b, p, lo, hi)
  m = size(A, 1)
  n = size(A, 2)
  A = full(A)

  nlo = 0
  for i = 1:length(lo)
    if lo[i] > 0
      nlo += 1
    end
  end
  nhi = 0
  for i = 1:length(hi)
    if hi[i] < 1e300
      nhi += 1
    end
  end

  A = [A; zeros(nlo + nhi, n)]
  b = [b; zeros(nlo + nhi)]

  nlo = 0
  for i = 1:length(lo)
    if lo[i] > 0
      nlo += 1
      A[m + nlo, i] = 1
      b[m + nlo] = lo[i]
    end
  end

  nhi = 0
  for i = 1:length(hi)
    if hi[i] < 1e300
      nhi += 1
      A[m + nlo + nhi, i] = 1
      b[m + nlo + nhi] = hi[i]
    end
  end

  new_n = length(b)
  A = [A [zeros(m, nlo + nhi); -eye(nlo) zeros(nlo, nhi); zeros(nhi, nlo) eye(nhi)]]
  A = SparseMatrixCSC{Float64, Int32}(sparse(A))
  p = [p; zeros(nlo + nhi)]

  A, b, p
end

function load_ipm_input(name)
  A = full(matrix_market_read(string(name, ".mtx")))'
  b = matrix_market_read(string(name, "_b.mtx"))
  p = matrix_market_read(string(name, "_c.mtx"))

  n = string(name, "_lo.mtx")
  if isfile(n)
    lo = matrix_market_read(string(name, "_lo.mtx"))
  else
    lo = zeros(0)
  end
  n = string(name, "_hi.mtx")
  if isfile(n)
    hi = matrix_market_read(string(name, "_hi.mtx"))
  else
    hi = zeros(0)
  end

  canonicalize_ipm_input(A, b, p, lo, hi)
end

type Factor <: Factorization{Float64}
  p::Vector{Int}
  function Factor(p::Vector{Int})
    new(p)
  end
end

function cholfact_int32(B :: SparseMatrixCSC{Float64, Int32})
    R = Factor(dss_analyze(B))
    dss_factor(R.p, B)
    R
end

import Base.LinAlg: (\)

function (\)(L::Factor, b::Vector)
    y = zeros(b)
    dss_solve!(L.p, b, y)
    y
end

# Compared to the Jongsoo's original ipm_ref(), the following changes are made
# to make it work around OptFramework issues:
#(1) No default parameter
#(2) No print
#(3) cholfact_int32(B) instead of cholfact(SparseMatrixCSC{Float64, Int64}(B))

# primal-dual interior-point method for problem
# min p'*x s.t. Ax = b, x >= 0

function ipm_ref(A, b, p) # A: constraint coefficients, b: constraint rhs, p: objective
  t0 = time()
  (m,n) = size(A)

# set initial point, based on largest element in (A,b,p)
  bigM = maximum(A)
  bigM = maximum([norm(b, Inf) norm(p, Inf) bigM])
  
  # Do not know why, Julia cannot automatically figure out the types of x, s and y,
  # for which we miss some patterns. Here manually annotate the types.
  x::Vector{Float64} = 100*bigM*ones(n)
  s::Vector{Float64} = x
  y::Vector{Float64} = zeros(m)

  bc = 1 + maximum([norm(b) norm(p)])

  relResidual = NaN
  iter=1

  blas1_time = 0.
  spgemm_time = 0.
  fact_time = 0.
  trslv_time = 0.
  spmv_time = 0.

  D = speye_int32(n)
  Rp = zeros(size(A, 1))
  
  # We need Rc to live into the loop so that we can create a temporary based on it.
  # TODO: remove this statement after a solution is found.
  Rc = zeros(x)

  for iter=1:200
    set_matrix_property(Dict(
        :Rc  => SA_HAS_DEDICATED_MEMORY,
      )
    )

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
    R = cholfact_int32(B)
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

  x, time() - t0, spgemm_time, fact_time, blas1_time, trslv_time, spmv_time, iter, relResidual, f[1][1]
end
