using SparseAccelerator

function abiatic(Hdmat, d, nqbits, T)
  total_time = -time()

  nsets = 2^nqbits
  Phi0at0 = 1/sqrt(2^nqbits)*ones(nsets)

  ode_time = -time()

  h = 0.01
  rtol = 1e-6

  # Compute the constants once
  c30 = 3/8
  c31 = 3/32
  c32 = 9/32
  c40 = 12/13
  c41 = 1932/2197
  c42 = -7200/2197
  c43 = 7296/2197
  c51 = 439/216
  c52 = -8
  c53 = 3680/513
  c54 = -845/4104
  c61 = -8/27
  c62 = 2
  c63 = -3544/2565
  c64 = 1859/4104
  c65 = -11/40
  cw1 = 25/216
  cw3 = 1408/2565
  cw4 = 2197/4104
  cw5 = -1/5
  cz1 = 16/135
  cz3 = 6656/12825
  cz4 = 28561/56430
  cz5 = -9/50
  cz6 = 2/55
  ce1 = 1/360
  ce3 = -128/4275
  ce4 = -2197/75240
  ce5 = 1/50
  ce6 = 2/55

  # Absolute tolerance
  atol = 1e-13
  alpha = 0.8
  k = 0
  # Initial time moment
  i = 1
  Tv = zeros(0)
  push!(Tv, 0)
  t = 0
  # Initial condition
  meanenergy = zeros(1)
  variance = zeros(1)
  meanenergy[1] = dot(vec(abs(Phi0at0).^2), d)
  variance[1] = dot(vec(abs(Phi0at0).^2), (d - meanenergy[1]).^2)
  wi = Phi0at0
  # If it is the last iteration, then lastit = 1, otherwise lastit = 0
  lastit = 0
  spmv_time = 0

  while lastit == 0
      # Stretch the step if within 10% of b-t
      if t + 1.1*h > T
          h = T - t
          lastit = 1
      end
     
      # Compute the step
      s = t/T
      w = wi
      spmv_time -= time()
      s1 = -im*(-(1 - s)*Hdmat*w + s*d.*w)
      spmv_time += time()

      s = (t + 0.25*h)/T
      w = wi + 0.25*h*s1
      spmv_time -= time()
      s2 = -im*(-(1 - s)*Hdmat*w + s*d.*w)
      spmv_time += time()

      s = (t + c30*h)/T
      w = wi + c31*h*s1 + c32*h*s2
      spmv_time -= time()
      s3 = -im*(-(1 - s)*Hdmat*w + s*d.*w)
      spmv_time += time()

      s = (t + c40*h)/T
      w = wi + c41*h*s1 + c42*h*s2 + c43*h*s3
      spmv_time -= time()
      s4 = -im*(-(1 - s)*Hdmat*w + s*d.*w)
      spmv_time += time()

      s = (t + h)/T
      w = wi + c51*h*s1 + c52*h*s2 + c53*h*s3 + c54*h*s4
      spmv_time -= time()
      s5 = -im*(-(1 - s)*Hdmat*w + s*d.*w)
      spmv_time += time()

      s = (t + 0.5*h)/T
      w = wi + c61*h*s1 + c62*h*s2 + c63*h*s3 + c64*h*s4 + c65*h*s5
      spmv_time -= time()
      s6 = -im*(-(1 - s)*Hdmat*w + s*d.*w)
      spmv_time += time()

      w = wi + h*(cw1*s1 + cw3*s3 + cw4*s4 + cw5*s5)
      z = wi + h*(cz1*s1 + cz3*s3 + cz4*s4 + cz5*s5 + cz6*s6)
      e = h * norm(ce1*s1 + ce3*s3 + ce4*s4 + ce5*s5 + ce6*s6)
      
      # Target tolerance for this step
      target_tol = rtol*norm(wi) + atol
      if e <= target_tol # In case the tolerance is met
          t = t + h
          h = alpha*h*(target_tol/e)^0.2
          i = i + 1
          push!(Tv, t)
          wi = z
          push!(meanenergy, dot(abs(z).^2, d))
          push!(variance, dot(abs(z).^2, (d - meanenergy[end]).^2))
          k = 0
      elseif k == 0 # Tolerance is not met for the first time in this step
          h = alpha*h*(target_tol/e)^0.2
          k = k + 1
          lastit = 0
      else # Tolerance is not met more than once in this step
          h = h/2
          lastit = 0
      end
  end

  println(size(Tv))

  ode_time += time()

  lanczos_time = -time()
  spectralgap = zeros(length(Tv))

  nit = min(30, nsets)

  x0 = rand(nsets)

  alpha = zeros(nit)
  beta = zeros(nit + 1)

  for i = 1:length(Tv)
    s = Tv[i]/T

    q0 = zeros(nsets)
    q1 = x0/norm(x0)

    for k = 1:nit
      spmv_time -= time()
      uk = -(1 - s)*Hdmat*q1 + s*d.*q1
      spmv_time += time()

      alpha[k] = dot(q1, uk)
      uk -= beta[k]*q0 + alpha[k]*q1
      beta[k + 1] = norm(uk)
      if beta[k + 1] == 0
        println("Error")
        return
      end

      q0 = q1
      q1 = uk/beta[k + 1]
    end

    TT = SymTridiagonal(alpha, beta[2:nit])

    Dt = sort(eig(TT)[1], 1)
    Dt = [Dt[1]; Dt[end]]

    spectralgap[i] = abs(Dt[1] - Dt[2])^2
  end
  lanczos_time += time()

  optimalenergy = minimum(d)

  total_time += time()

  @printf("%10s%10s%12s%12s%12s\n", "time", "eigengap", "optimal", "meanenergy", "variance");
  for i = length(Tv):length(Tv)
    @printf("%10.2f%10.4f%12.4f%12.4f%12.4f\n", Tv[i], spectralgap[i], optimalenergy, meanenergy[i], variance[i])
  end

  println("total_time = $(total_time) ode_time = $(ode_time) lanczos_time = $(lanczos_time)")
  println("spmv_time = $spmv_time")
end

function spmv_test(A)
  mknobA = (SparseAccelerator.new_matrix_knob)(A, true, true, true, true, false, false)
  fknob_spmv = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobA, fknob_spmv)

  x = rand(size(A, 1))
  y = zeros(size(A, 1))

  for i = 1:100
    SparseAccelerator.SpMV!(y, A, x, fknob_spmv)

    temp = y
    y = x
    x = temp
  end
  
  x
end

function spmv_test_reordering(A)
  mknobA = (SparseAccelerator.new_matrix_knob)(A, true, true, true, true, false, false)
  fknob_spmv = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobA, fknob_spmv)

  (SparseAccelerator.set_reordering_decision_maker)(fknob_spmv)

  x = rand(size(A, 1))
  y = zeros(size(A, 1))

  for i = 1:100
    SparseAccelerator.SpMV!(y, A, x, fknob_spmv)

    temp = y
    y = x
    x = temp
  end
  
  x
end

nqbits = parse(Int, ARGS[1])
T = parse(Int, ARGS[2])

nsets = 2^nqbits
d = 4*collect(-nsets/2 : nsets/2 - 1)

tunnelstr = -1

hdmat_time = -time()
Hdmat = zeros(nsets, nsets)
for i=1:nsets
  for j=0:nqbits-1
    flip = (i - 1)$(1 << j) # flip jth bit in zero-based indexing
    Hdmat[i,flip+1] = 1
  end
end
Hdmat = SparseMatrixCSC{Float64, Int32}(sparse(Hdmat))

#ccall((:store_matrix_market, SparseAccelerator.LIB_PATH), Void,
       #(Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
       #ARGS[3], Hdmat.m, Hdmat.n, Hdmat.colptr, Hdmat.rowval, Hdmat.nzval)

abiatic(Hdmat, d, nqbits, T)
