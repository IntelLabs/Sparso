include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)

function abiatic(Hdmat, d, nqbits, T)
  total_time = -time()

  nsets = 2^nqbits
  Phi0at0 = 1/sqrt(2^nqbits)*ones(nsets)

  ode_time = -time()

  h = 0.01
  rtol = 1e-6

  # Compute the constants once
  c30 = 3/8
  c31 = -im*3/32
  c32 = -im*9/32
  c40 = 12/13
  c41 = -im*1932/2197
  c42 = im*7200/2197
  c43 = -im*7296/2197
  c51 = -im*439/216
  c52 = im*8
  c53 = -im*3680/513
  c54 = im*845/4104
  c61 = im*8/27
  c62 = -im*2
  c63 = im*3544/2565
  c64 = -im*1859/4104
  c65 = im*11/40
  cz1 = -im*16/135
  cz3 = -im*6656/12825
  cz4 = -im*28561/56430
  cz5 = im*9/50
  cz6 = -im*2/55
  ce1 = -im*1/360
  ce3 = im*128/4275
  ce4 = im*2197/75240
  ce5 = -im*1/50
  ce6 = -im*2/55

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
      s = t/T # real
      w = wi
      spmv_time -= time()
      s1 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + 0.25*h)/T
      w = wi + -im*0.25*h*s1
      spmv_time -= time()
      s2 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + c30*h)/T
      w = wi + c31*h*s1 + c32*h*s2
      spmv_time -= time()
      s3 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + c40*h)/T
      w = wi + c41*h*s1 + c42*h*s2 + c43*h*s3
      spmv_time -= time()
      s4 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + h)/T
      w = wi + c51*h*s1 + c52*h*s2 + c53*h*s3 + c54*h*s4
      spmv_time -= time()
      s5 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + 0.5*h)/T
      w = wi + c61*h*s1 + c62*h*s2 + c63*h*s3 + c64*h*s4 + c65*h*s5
      spmv_time -= time()
      s6 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      z = wi + h*(cz1*s1 + cz3*s3 + cz4*s4 + cz5*s5 + cz6*s6)
      e = h * norm(ce1*s1 + ce3*s3 + ce4*s4 + ce5*s5 + ce6*s6) # real
      
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

  println("total_time = $(total_time) ode_time = $(ode_time) lanczos_time = $(lanczos_time) spmv_time = $(spmv_time)")
end

function print_result(Tv, spectralgap, optimalenergy, meanenergy, variance)
  @printf("%10s%10s%12s%12s%12s\n", "time", "eigengap", "optimal", "meanenergy", "variance");
  for i = length(Tv):length(Tv)
    @printf("%10.2f%10.4f%12.4f%12.4f%12.4f\n", Tv[i], spectralgap[i], optimalenergy, meanenergy[i], variance[i])
  end
end

function abiatic_sa(Hdmat, d, nqbits, T)
  total_time = -time()

  set_matrix_property(Dict(
     :Hdmat => SA_SYMM_VALUED | SA_STRUCTURE_ONLY, 
     )
  )

  nsets = 2^nqbits
  Phi0at0 = 1/sqrt(2^nqbits)*ones(nsets)

  ode_time = -time()

  h = 0.01
  rtol = 1e-6

  # Compute the constants once
  c30 = 3/8
  c31 = -im*3/32
  c32 = -im*9/32
  c40 = 12/13
  c41 = -im*1932/2197
  c42 = im*7200/2197
  c43 = -im*7296/2197
  c51 = -im*439/216
  c52 = im*8
  c53 = -im*3680/513
  c54 = im*845/4104
  c61 = im*8/27
  c62 = -im*2
  c63 = im*3544/2565
  c64 = -im*1859/4104
  c65 = im*11/40
  cz1 = -im*16/135
  cz3 = -im*6656/12825
  cz4 = -im*28561/56430
  cz5 = im*9/50
  cz6 = -im*2/55
  ce1 = -im*1/360
  ce3 = im*128/4275
  ce4 = im*2197/75240
  ce5 = -im*1/50
  ce6 = -im*2/55

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
  wi = Array{Complex128}(Phi0at0)
  # If it is the last iteration, then lastit = 1, otherwise lastit = 0
  lastit = 0
  spmv_time = 0

  s1 = Array{Complex128}(length(d))
  s2 = Array{Complex128}(length(d))
  s3 = Array{Complex128}(length(d))
  s4 = Array{Complex128}(length(d))
  s5 = Array{Complex128}(length(d))
  s6 = Array{Complex128}(length(d))

  while lastit == 0
      # Stretch the step if within 10% of b-t
      if t + 1.1*h > T
          h = T - t
          lastit = 1
      end
     
      # Compute the step
      s = t/T # real
      w = wi
      spmv_time -= time()
      s1 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + 0.25*h)/T
      w = wi + -im*0.25*h*s1
      spmv_time -= time()
      s2 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + c30*h)/T
      w = wi + c31*h*s1 + c32*h*s2
      spmv_time -= time()
      s3 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + c40*h)/T
      w = wi + c41*h*s1 + c42*h*s2 + c43*h*s3
      spmv_time -= time()
      s4 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + h)/T
      w = wi + c51*h*s1 + c52*h*s2 + c53*h*s3 + c54*h*s4
      spmv_time -= time()
      s5 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      s = (t + 0.5*h)/T
      w = wi + c61*h*s1 + c62*h*s2 + c63*h*s3 + c64*h*s4 + c65*h*s5
      spmv_time -= time()
      s6 = -(1 - s)*Hdmat*w + s*d.*w
      spmv_time += time()

      z = wi + h*(cz1*s1 + cz3*s3 + cz4*s4 + cz5*s5 + cz6*s6)
      e = h * norm(ce1*s1 + ce3*s3 + ce4*s4 + ce5*s5 + ce6*s6) # real
      
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

  uk = Array{Float64}(length(d))

  dot_time = 0
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
      #q1 = uk/beta[k + 1]
      q1 = SparseAccelerator.WAXPB(1/beta[k+1], uk, 0) 
    end

    TT = SymTridiagonal(alpha, beta[2:nit])

    Dt = sort(eig(TT)[1], 1)
    Dt = [Dt[1]; Dt[end]]

    spectralgap[i] = abs(Dt[1] - Dt[2])^2
  end
  lanczos_time += time()

  optimalenergy = minimum(d)

  total_time += time()

  print_result(Tv, spectralgap, optimalenergy, meanenergy, variance)

  println("total_time = $(total_time) ode_time = $(ode_time) lanczos_time = $(lanczos_time) spmv_time = $(spmv_time)")
end


function abiatic_opt(Hdmat, d, nqbits, T, with_context_opt)
  total_time = -time()

  if with_context_opt
    __mknobHdmat = (SparseAccelerator.new_matrix_knob)(Hdmat, true, true, true, true, true, false)

    fknob_spmv1 = (SparseAccelerator.new_function_knob)()
    fknob_spmv2 = (SparseAccelerator.new_function_knob)()
    fknob_spmv3 = (SparseAccelerator.new_function_knob)()
    fknob_spmv4 = (SparseAccelerator.new_function_knob)()
    fknob_spmv5 = (SparseAccelerator.new_function_knob)()
    fknob_spmv6 = (SparseAccelerator.new_function_knob)()
    fknob_spmv7 = (SparseAccelerator.new_function_knob)()

    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv1)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv2)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv3)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv4)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv5)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv6)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv7)
  else
    fknob_spmv1 = C_NULL 
    fknob_spmv2 = C_NULL
    fknob_spmv3 = C_NULL
    fknob_spmv4 = C_NULL
    fknob_spmv5 = C_NULL
    fknob_spmv6 = C_NULL
    fknob_spmv7 = C_NULL
  end

  nsets = 2^nqbits
  Phi0at0 = 1/sqrt(2^nqbits)*ones(nsets)

  ode_time = -time()

  h = 0.01
  rtol = 1e-6

  # Compute the constants once
  c30 = 3/8
  c31 = -im*3/32
  c32 = -im*9/32
  c40 = 12/13
  c41 = -im*1932/2197
  c42 = im*7200/2197
  c43 = -im*7296/2197
  c51 = -im*439/216
  c52 = im*8
  c53 = -im*3680/513
  c54 = im*845/4104
  c61 = im*8/27
  c62 = -im*2
  c63 = im*3544/2565
  c64 = -im*1859/4104
  c65 = im*11/40
  cz1 = -im*16/135
  cz3 = -im*6656/12825
  cz4 = -im*28561/56430
  cz5 = im*9/50
  cz6 = -im*2/55
  ce1 = -im*1/360
  ce3 = im*128/4275
  ce4 = im*2197/75240
  ce5 = -im*1/50
  ce6 = -im*2/55

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
  wi = Array{Complex128}(Phi0at0)
  # If it is the last iteration, then lastit = 1, otherwise lastit = 0
  lastit = 0
  spmv_time = 0

  w = Array{Complex128}(length(d))
  temp = Array{Complex128}(length(d))
  s1 = Array{Complex128}(length(d))
  s2 = Array{Complex128}(length(d))
  s3 = Array{Complex128}(length(d))
  s4 = Array{Complex128}(length(d))
  s5 = Array{Complex128}(length(d))
  s6 = Array{Complex128}(length(d))
  z = Array{Complex128}(length(d))

  while lastit == 0
      # Stretch the step if within 10% of b-t
      if t + 1.1*h > T
          h = T - t
          lastit = 1
      end
     
      # Compute the step
      s = t/T # real
      #w = copy(wi)
      SparseAccelerator.copy!(w, wi)
      spmv_time -= time()
      #s1 = -(1 - s)*Hdmat*w + s*d.*w
      SparseAccelerator.element_wise_multiply!(temp, d, w)
      SparseAccelerator.SpMV!(s1, s - 1, Hdmat, w, s, temp, fknob_spmv1)
      spmv_time += time()

      s = (t + 0.25*h)/T
      #w = wi + -im*0.25*h*s1
      SparseAccelerator.WAXPBY!(w, 1, wi, -im*0.25*h, s1)
      spmv_time -= time()
      #s2 = -(1 - s)*Hdmat*w + s*d.*w
      SparseAccelerator.element_wise_multiply!(temp, d, w)
      SparseAccelerator.SpMV!(s2, s - 1, Hdmat, w, s, temp, fknob_spmv2)
      spmv_time += time()

      s = (t + c30*h)/T
      #w = wi + c31*h*s1 + c32*h*s2
      SparseAccelerator.WAXPBY!(w, 1, wi, c31*h, s1)
      SparseAccelerator.WAXPBY!(w, 1, w, c32*h, s2)
      spmv_time -= time()
      #s3 = -(1 - s)*Hdmat*w + s*d.*w
      SparseAccelerator.element_wise_multiply!(temp, d, w)
      SparseAccelerator.SpMV!(s3, s - 1, Hdmat, w, s, temp, fknob_spmv3)
      spmv_time += time()

      s = (t + c40*h)/T
      #w = wi + c41*h*s1 + c42*h*s2 + c43*h*s3
      SparseAccelerator.WAXPBY!(w, 1, wi, c41*h, s1)
      SparseAccelerator.WAXPBY!(w, 1, w, c42*h, s2)
      SparseAccelerator.WAXPBY!(w, 1, w, c43*h, s3)
      spmv_time -= time()
      #s4 = -(1 - s)*Hdmat*w + s*d.*w
      SparseAccelerator.element_wise_multiply!(temp, d, w)
      SparseAccelerator.SpMV!(s4, s - 1, Hdmat, w, s, temp, fknob_spmv4)
      spmv_time += time()

      s = (t + h)/T
      #w = wi + c51*h*s1 + c52*h*s2 + c53*h*s3 + c54*h*s4
      SparseAccelerator.WAXPBY!(w, 1, wi, c51*h, s1)
      SparseAccelerator.WAXPBY!(w, 1, w, c52*h, s2)
      SparseAccelerator.WAXPBY!(w, 1, w, c53*h, s3)
      SparseAccelerator.WAXPBY!(w, 1, w, c54*h, s4)
      spmv_time -= time()
      #s5 = -(1 - s)*Hdmat*w + s*d.*w
      SparseAccelerator.element_wise_multiply!(temp, d, w)
      SparseAccelerator.SpMV!(s5, s - 1, Hdmat, w, s, temp, fknob_spmv5)
      spmv_time += time()

      s = (t + 0.5*h)/T
      #w = wi + c61*h*s1 + c62*h*s2 + c63*h*s3 + c64*h*s4 + c65*h*s5
      SparseAccelerator.WAXPBY!(w, 1, wi, c61*h, s1)
      SparseAccelerator.WAXPBY!(w, 1, w, c62*h, s2)
      SparseAccelerator.WAXPBY!(w, 1, w, c63*h, s3)
      SparseAccelerator.WAXPBY!(w, 1, w, c64*h, s4)
      SparseAccelerator.WAXPBY!(w, 1, w, c65*h, s5)
      spmv_time -= time()
      #s6 = -(1 - s)*Hdmat*w + s*d.*w
      SparseAccelerator.element_wise_multiply!(temp, d, w)
      SparseAccelerator.SpMV!(s6, s - 1, Hdmat, w, s, temp, fknob_spmv6)
      spmv_time += time()

      #z = wi + h*(cz1*s1 + cz3*s3 + cz4*s4 + cz5*s5 + cz6*s6)
      SparseAccelerator.WAXPBY!(w, cz1, s1, cz3, s3)
      SparseAccelerator.WAXPBY!(w, 1, w, cz4, s4)
      SparseAccelerator.WAXPBY!(w, 1, w, cz5, s5)
      SparseAccelerator.WAXPBY!(w, 1, w, cz6, s6)
      SparseAccelerator.WAXPBY!(z, 1, wi, h, w)
      #e = h * norm(ce1*s1 + ce3*s3 + ce4*s4 + ce5*s5 + ce6*s6) # real
      SparseAccelerator.WAXPBY!(w, ce1, s1, ce3, s3)
      SparseAccelerator.WAXPBY!(w, 1, w, ce4, s4)
      SparseAccelerator.WAXPBY!(w, 1, w, ce5, s5)
      SparseAccelerator.WAXPBY!(w, 1, w, ce6, s6)
      e = h*SparseAccelerator.norm(w)
      
      # Target tolerance for this step
      target_tol = rtol*SparseAccelerator.norm(wi) + atol
      if e <= target_tol # In case the tolerance is met
          t = t + h
          h = alpha*h*(target_tol/e)^0.2
          i = i + 1
          push!(Tv, t)
          SparseAccelerator.copy!(wi, z)
          push!(meanenergy, SparseAccelerator.dot(SparseAccelerator.abs(z).^2, d))
          push!(variance, SparseAccelerator.dot(SparseAccelerator.abs(z).^2, (d - meanenergy[end]).^2))
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

  uk = Array{Float64}(nsets)

  temp2 = Array{Float64}(nsets)
  temp3 = Array{Float64}(nsets)
  temp4 = Array{Float64}(nsets)
  temp5 = Array{Float64}(nsets)

  for i = 1:length(Tv)
    s = Tv[i]/T

    q0 = zeros(nsets)
    #q1 = x0/SparseAccelerator.norm(x0)
    SparseAccelerator.WAXPB!(temp5, 1/SparseAccelerator.norm(x0), x0, 0) 
    q1 = temp5

    for k = 1:nit
      spmv_time -= time()
      #uk = -(1 - s)*Hdmat*q1 + s*d.*q1
      SparseAccelerator.element_wise_multiply!(temp2, d, q1)
      SparseAccelerator.SpMV!(uk, s - 1, Hdmat, q1, s, temp2, 0, fknob_spmv7)
      spmv_time += time()

      #alpha[k] = dot(q1, uk)
      alpha[k] = SparseAccelerator.dot(q1, uk)

      #uk -= beta[k]*q0 + alpha[k]*q1
      #SparseAccelerator.WAXPBY!(uk, 1, uk, -beta[k], q0)
      #SparseAccelerator.WAXPBY!(uk, 1, uk, -alpha[k], q1)
      SparseAccelerator.WAXPBY!(temp3, beta[k], q0, alpha[k], q1)
      SparseAccelerator.WAXPBY!(uk, 1, uk, -1, temp3)

      #beta[k + 1] = norm(uk)
      beta[k + 1] = SparseAccelerator.norm(uk)

      if beta[k + 1] == 0
        println("Error")
        return
      end

      q0 = q1
      #q1 = uk/beta[k + 1]
      SparseAccelerator.WAXPB!(temp5, 1/beta[k+1], uk, 0) 
      q1 = temp5
    end

    TT = SymTridiagonal(alpha, beta[2:nit])

    Dt = sort(eig(TT)[1], 1)
    Dt = [Dt[1]; Dt[end]]

    spectralgap[i] = abs(Dt[1] - Dt[2])^2
  end
  lanczos_time += time()

  optimalenergy = minimum(d)

  if with_context_opt
    (SparseAccelerator.delete_matrix_knob)(__mknobHdmat)

    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv1)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv2)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv3)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv4)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv5)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv6)
    (SparseAccelerator.add_mknob_to_fknob)(__mknobHdmat, fknob_spmv7)
  end

  total_time += time()

  @printf("%10s%10s%12s%12s%12s\n", "time", "eigengap", "optimal", "meanenergy", "variance");
  for i = length(Tv):length(Tv)
    @printf("%10.2f%10.4f%12.4f%12.4f%12.4f\n", Tv[i], spectralgap[i], optimalenergy, meanenergy[i], variance[i])
  end

  println("total_time = $(total_time) ode_time = $(ode_time) lanczos_time = $(lanczos_time) spmv_time = $(spmv_time)")
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

if length(ARGS) == 3
  tests = ARGS[3]
else
  tests = "OMS"
end

if contains(tests, "O")
  println("\nOriginal: ")
  srand(0)
  println(@time(abiatic(Hdmat, d, nqbits, T)))
  println("\nEND Oringinal.")
end

#SparseAccelerator.set_knob_log_level(1)

if contains(tests, "M")
  println("\nManual_call_replacement:")
  srand(0)
  println(@time(abiatic_opt(Hdmat, d, nqbits, T, false)))
  println("\nEND Manual_call_replacement.")

  println("\nManual_call_replacement_and_context_opt:")
  srand(0)
  println(@time(abiatic_opt(Hdmat, d, nqbits, T, true)))
  println("\nEND Manual_call_replacement_and_context_opt.")
end


if contains(tests, "S")
#println("\nSparso_call_replacement:")
#srand(0)
#println(@time(@acc abiatic_sa(Hdmat, d, nqbits, T, false)))
#println("\nEND Sparso_call_replacement.")

  println("\nSparso_call_replacement_and_context_opt:")
  srand(0)
  println(@time(@acc abiatic_sa(Hdmat, d, nqbits, T)))
  println("\nEND Sparso_call_replacement_and_context_opt.")
end
