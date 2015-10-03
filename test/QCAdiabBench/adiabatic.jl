using SparseAccelerator

# Runge-Kutta-Fehberg order 4/order 5 embedded pair
# Input:
# f - inline function f(t, y)
# a, b - interval
# ya - initial condition
# h - initial step-size
# rtol - relative error tolerance per step
#
# Output:
# w - computed solution
# t - time steps
#
# Examples:
# [y t]=rkf45(@myfunc,0,1,1,0.1,1e-4);          here 'myfunc' is a user-defined function in M-file
# [w t]=rkf45(inline('sin(y*t)','t','y'),0,1,1,0.1,0.01);
# f=inline('sin(y(1))-cos(y(2))','t','y');
# [y t]=rkf45(f,0,1,1,0.01,1e-5);

function rkf45(a, b, ya, h, rtol, Hdmat, d)
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
  tt = zeros(0)
  push!(tt, a)
  t = a
  # Initial condition
  y = ya'
  wi = ya
  # If it is the last iteration, then lastit = 1, otherwise lastit = 0
  lastit = 0
  time_spmv = 0

  while lastit == 0
      # Stretch the step if within 10% of b-t
      if t + 1.1*h > b
          h = b - t
          lastit = 1
      end
     
      # Compute the step
      s = t/b
      w = wi
      time_spmv -= time()
      s1 = -im*((1 - s)*Hdmat*w + s*d.*w)
      time_spmv += time()

      s = (t + 0.25*h)/b
      w = wi + 0.25*h*s1
      time_spmv -= time()
      s2 = -im*((1 - s)*Hdmat*w + s*d.*w)
      time_spmv += time()

      s = (t + c30*h)/b
      w = wi + c31*h*s1 + c32*h*s2
      time_spmv -= time()
      s3 = -im*((1 - s)*Hdmat*w + s*d.*w)
      time_spmv += time()

      s = (t + c40*h)/b
      w = wi + c41*h*s1 + c42*h*s2 + c43*h*s3
      time_spmv -= time()
      s4 = -im*((1 - s)*Hdmat*w + s*d.*w)
      time_spmv += time()

      s = (t + h)/b
      w = wi + c51*h*s1 + c52*h*s2 + c53*h*s3 + c54*h*s4
      time_spmv -= time()
      s5 = -im*((1 - s)*Hdmat*w + s*d.*w)
      time_spmv += time()

      s = (t + 0.5*h)/b
      w = wi + c61*h*s1 + c62*h*s2 + c63*h*s3 + c64*h*s4 + c65*h*s5
      time_spmv -= time()
      s6 = -im*((1 - s)*Hdmat*w + s*d.*w)
      time_spmv += time()

      w = wi + h * (cw1 * s1 + cw3 * s3 + cw4 * s4 + cw5 * s5)
      z = wi + h * (cz1 * s1 + cz3 * s3 + cz4 * s4 + cz5 * s5 + cz6 * s6)
      e = h * norm(ce1 * s1 + ce3 * s3 + ce4 * s4 + ce5 * s5 + ce6 * s6)
      
      # Target tolerance for this step
      T = rtol*norm(wi) + atol
      if e <= T # In case the tolerance is met
          t = t + h
          h = alpha*h*(T/e)^0.2
          i = i + 1
          push!(tt, t)
          wi = z
          y = [y; z']
          k = 0
      elseif k == 0 # Tolerance is not met for the first time in this step
          h = alpha*h*(T/e)^0.2
          k = k + 1
          lastit = 0
      else # Tolerance is not met more than once in this step
          h = h/2
          lastit = 0
      end
  end

  println(size(tt))
  println("time_spmv = $time_spmv")

  y, tt
end

function lanczos(A, K, nit, check)
  d = size(A)
  assert(nit <= d[1])
  assert(rem(K,2) == 0) # simplifes logic a bit

  n = d[1]
  q0 = zeros(n,1)
  alpha = zeros(0)
  beta = zeros(1)
  x0 = rand(n,1)
  #println(x0)
  q1 = x0/norm(x0)
  Q = [q0 q1]
  #println("A = $A")
  #println("Q = $Q")
  for k = 1:nit
      uk = A*Q[:, k + 1]
      #println("uk = $uk")
      push!(alpha, dot(Q[:, k+1], uk))
      uk -= beta[k]*Q[:, k] + alpha[k]*Q[:, k+1]
      #println("uk = $uk")
      push!(beta, norm(uk))
      if beta[k + 1] == 0
          println("Error")
          return
      end
      Q = [Q uk/beta[k+1]]
  end

  #println(alpha)
  #println(beta)
  T = diagm(alpha) + diagm(beta[2:nit], -1) + diagm(beta[2:nit], 1);
  #println("T = $T")

  Dt = sort(eig(T)[1], 1)
  Ko2 = Int32(round(K/2))
  Dt = [Dt[1:Ko2]; Dt[end - Ko2 + 1:end]]
  if check
    E = sort(eigs(A, K, which="BE")', 2)
    error = abs(Dt-E)./abs(E)
    println(error)
  end

  Dt
end

function abiatic(d, nqbits, T)
  total_time = -time()

  nsets = 2^nqbits
  tunnelstr = 1

  hdmat_time = -time()
  Hdmat = zeros(nsets, nsets)
  for i=1:nsets
    for j=0:nqbits-1
      flip = (i - 1)$(1 << j) # flip jth bit in zero-based indexing
      Hdmat[i,flip+1] = -tunnelstr
    end
  end
  Hdmat = sparse(Hdmat)

  Hpmat = spdiagm(d)
  Phi0at0 = 1/sqrt(2^nqbits)*ones(nsets)
  hdmat_time += time()

  ode_time = -time()
  Y, Tv = rkf45(0, T, Phi0at0, 0.01, 1e-6, Hdmat, d)
  ode_time += time()

  lanczos_time = -time()
  spectralgap = zeros(0)
  for i = 1:length(Tv)
    s = Tv[i]/T
    H = (1-s)*Hdmat + s*Hpmat
    lambdat = lanczos(H, 2, min(30, nsets), false)
    push!(spectralgap, abs(lambdat[1] - lambdat[2])^2)
  end
  lanczos_time += time()

  result_time = -time()
  meanenergy = abs(Y).^2*diag(Hpmat)
  variance = zeros(length(meanenergy))
  for i=1:length(meanenergy)
    variance[i] = dot(vec(abs(Y[i, :]).^2), (diag(Hpmat) - meanenergy[i]).^2)
  end
  optimalenergy = minimum(Hpmat)
  result_time += time()

  total_time += time()

  @printf("%10s%10s%12s%12s%12s\n", "time", "eigengap", "optimal", "meanenergy", "variance");
  for i = length(Tv):length(Tv)
    @printf("%10.2f%10.4f%12.4f%12.4f%12.4f\n", Tv[i], spectralgap[i], optimalenergy, meanenergy[i], variance[i])
  end

  println("total_time = $(total_time) hdmat_time = $(hdmat_time) ode_time = $(ode_time) lanczos_time = $(lanczos_time) result_time = $(result_time)")
end

nqbits = parse(Int, ARGS[1])
T = parse(Int, ARGS[2])

nsets = 2^nqbits
d = 4*collect(-nsets/2 : nsets/2 - 1)

abiatic(d, nqbits, T)
