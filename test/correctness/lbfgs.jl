include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)

function mylogsumexp(b)
  # does logsumexp across column
  log(1 + exp(-abs(b)))+max(b,0)
end

function LogisticLoss(w, X, Xt, y, lambda)
  Xw = X*w
  yXw = y.*Xw
  m = size(X, 1)

  s = sum(mylogsumexp(-yXw))
  fk = s/m + (lambda/2)*norm(w)^2
  gk = -(Xt*(y./(1+exp(yXw))))/m + lambda*w

  fk, gk
end

function GD(X, Xt, y, lambda, xinit, tol)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  t0 = time()
  spmv_time = 0

  x = xinit

  it = 1
  for it=1:100
    (n,p) = size(X)
    tic()
    Xw = X*x
    spmv_time += toq()
    yXw = y.*Xw

    s = sum(mylogsumexp(-yXw))
    fk0 = s/n+(lambda/2)*norm(x)^2
    dfk = -(Xt*(y./(1+exp(yXw))))/n + lambda*x

    if (norm(dfk) < tol)
      break;
    end

    # backtracking line search using armijo criterion
    alphaMax = 1 # this is the maximum step length
    alpha = alphaMax
    rho = 1/2 # < 1 reduction factor of alpha
    c_1 = 1e-4

    while true
      # logistic loss objective funcion
      w = x - alpha*dfk

      tic()
      Xw = X*w
      spmv_time += toq()
      yXw = y.*Xw

      s = sum(mylogsumexp(-yXw))
      fk = s/n+(lambda/2)*norm(w)^2
      # end objective function

      if (fk <= fk0 - c_1*alpha*norm(dfk)^2)
        break
      end

      alpha = rho*alpha

      if alpha < 10*eps()
        alpha=0.1
        error("Error in Line search - alpha close to working precision")
      end
    end
    # end backtracking line search using armijo criterion

    #@printf("[%5d]%10.3e%10.3e%10.7f\n", it, alpha, norm(dfk), fk0)
    x=x-alpha*dfk
  end
  println("SpMV takes $spmv_time sec.")
  println("GD takes $(time() - t0) sec.")
  x, it
end

function lbfgs_ref(X, y, lambda, xinit, tol, k)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  spmv_time = 0

  m, n = size(X)
  x = xinit
  Xt = X'

  S = zeros(n, k)
  Y = zeros(n, k)

  a = zeros(k, 1)

  t0 = time()
  it = 1

  # Pre-allocation of space. Should get rid of once Linxiang's object removal works
  # TODO: remove it.
#  Xw    = Array(Cdouble, length(y))
#  yXw   = Array(Cdouble, length(y))
#  dfkp1 = Array(Cdouble, length(x)) 
  temp  = Array(Cdouble, m)
  Xw    = Array(Cdouble, m)
  yXw   = Array(Cdouble, m)
  dfk   = Array(Cdouble, n)
  dfkp1 = Array(Cdouble, n) 
  w     = Array(Cdouble, n)
  dx    = Array(Cdouble, n)
  
#  temp = zeros(m)
#  Xw = zeros(m)
#  yXw = zeros(m)
#  dfk = zeros(n)
#  w = zeros(n)
#  dx = zeros(n)

  for it=1:100
    set_matrix_property(:Xt, SA_TRANSPOSE_OF, :X) 
    set_matrix_property(Dict(
        :temp  => SA_HAS_DEDICATED_MEMORY,
        :Xw    => SA_HAS_DEDICATED_MEMORY,
        :yXw   => SA_HAS_DEDICATED_MEMORY,
        :dfk   => SA_HAS_DEDICATED_MEMORY,
        :dfkp1 => SA_HAS_DEDICATED_MEMORY,
        :w     => SA_HAS_DEDICATED_MEMORY,
        :dx    => SA_HAS_DEDICATED_MEMORY,
      )
    )

    spmv_time -= time()
    #Xw = X*x
    A_mul_B!(Xw, X, x)
    spmv_time += time()

    yXw = y.*Xw

    s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
    fk0 = s/m+(lambda/2)*norm(x)^2
    temp = y./(1+exp(yXw))

    spmv_time -= time()
    dfk = -(Xt*temp)/m + lambda*x
    spmv_time += time()

    if (norm(dfk) < tol)
      break;
    end

    # Since Julia code for this is too slow we do this in C to show sufficient speedups by
    # faster SpMV. To be fair, the reference code also uses it.
    SparseAccelerator.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)

    # backtracking line search using armijo criterion
    alphaMax = 1 # this is the maximum step length
    alpha = alphaMax
    rho = 1/2 # < 1 reduction factor of alpha
    c_1 = 1e-4

    while true
      # logistic loss objective funcion
      w = x - alpha*dfk

      spmv_time -= time()
      Xw = X*w
      spmv_time += time()
      yXw = y.*Xw

      s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
      fk = s/m+(lambda/2)*norm(w)^2
      # end objective function

      if (fk <= fk0 - c_1*alpha*norm(dfk)^2)
        break
      end

      alpha = rho*alpha

      if alpha < 10*eps()
        alpha=0.1
        error("Error in Line search - alpha close to working precision")
      end
    end
    # end backtracking line search using armijo criterion

    #@printf("[%5d]%10.3e%10.3e%10.7f\n", it, alpha, norm(dfk), fk0)

    x = x + alpha*dx

    spmv_time -= time()
    Xw = X*x
    spmv_time += time()

    yXw = y.*Xw

    temp = y./(1+exp(yXw))

    spmv_time -= time()
    dfkp1 = -(Xt*temp)/m + lambda*x
    spmv_time += time()

    S[:,(it - 1)%k + 1] = alpha*dx
    Y[:,(it - 1)%k + 1] = dfkp1 - dfk
  end

  println("\nSpMV takes $spmv_time sec.")
  println("lbfgs_ref takes $(time() - t0) sec.")
  x, it
end

function lbfgs_opt(X, y, lambda, xinit, tol, k)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  spmv_time = 0
  log_time = 0
  direction_time = 0

  m, n = size(X)
  x = xinit
  Xt = X'

  spmv_count = 0

  mknobX = (SparseAccelerator.new_matrix_knob)(X, true, true, false, false, false, false)
  mknobXT = (SparseAccelerator.new_matrix_knob)(Xt, true, true, false, false, false, false)

  (SparseAccelerator.set_derivative)(mknobX, SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE, mknobXT)

  fknob_spmv1 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spmv1)
  fknob_spmv2 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobXT, fknob_spmv2)
  fknob_spmv3 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spmv3)
  fknob_spmv4 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spmv4)
  fknob_spmv5 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobXT, fknob_spmv5)

  S = zeros(n, k)
  Y = zeros(n, k)

  a = zeros(k, 1)

  dx    = Array(Cdouble, n)

  t0 = time()
  it = 1
  for it=1:100
    spmv_time -= time()
    #Xw = X*x
    Xw = SparseAccelerator.SpMV(X, x, fknob_spmv1)
    spmv_count += 1
    spmv_time += time()

    #yXw = y.*Xw
    yXw = SparseAccelerator.element_wise_multiply(y, Xw)

    log_time -= time()
    #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
    temp = zeros(length(yXw))
    SparseAccelerator.abs!(temp, yXw)
    SparseAccelerator.WAXPBY!(temp, -1, temp, 0, temp)
    SparseAccelerator.exp!(temp, temp)
    SparseAccelerator.log1p!(temp, temp)
    temp2 = zeros(length(yXw))
    SparseAccelerator.min!(temp2, yXw, 0)
    SparseAccelerator.WAXPBY!(temp, 1, temp, 1, temp2)
    s = SparseAccelerator.sum(temp)

    fk0 = s/m+(lambda/2)*SparseAccelerator.dot(x, x)
    #temp = y./(1+exp(yXw))
    SparseAccelerator.exp!(temp, yXw)
    SparseAccelerator.WAXPB!(temp, 1, temp, 1)
    SparseAccelerator.element_wise_divide!(temp, y, temp)

    log_time += time()

    spmv_time -= time()
    #dfk = -(Xt*temp)/m + lambda*x
    dfk = SparseAccelerator.SpMV(-1/m, Xt, temp, lambda, x, fknob_spmv2)
    spmv_count += 1
    spmv_time += time()

    if (SparseAccelerator.norm(dfk) < tol)
      break;
    end

    direction_time -= time()
    SparseAccelerator.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)
    direction_time += time()

    # backtracking line search using armijo criterion
    alphaMax = 1 # this is the maximum step length
    alpha = alphaMax
    rho = 1/2 # < 1 reduction factor of alpha
    c_1 = 1e-4

    while true
      # logistic loss objective funcion
      #w = x - alpha*dfk
      w = SparseAccelerator.WAXPBY(1, x, -alpha, dfk)

      spmv_time -= time()
      #Xw = X*w
      Xw = SparseAccelerator.SpMV(X, w, fknob_spmv3)
      spmv_count += 1
      spmv_time += time()
      SparseAccelerator.element_wise_multiply!(yXw, y, Xw)

      log_time -= time()
      #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
      SparseAccelerator.abs!(temp, yXw)
      SparseAccelerator.WAXPBY!(temp, -1, temp, 0, temp)
      SparseAccelerator.exp!(temp, temp)
      SparseAccelerator.log1p!(temp, temp)
      temp2 = zeros(length(yXw))
      SparseAccelerator.min!(temp2, yXw, 0)
      SparseAccelerator.WAXPBY!(temp, 1, temp, 1, temp2)
      s = SparseAccelerator.sum(temp)

      fk = s/m+(lambda/2)*SparseAccelerator.dot(w, w)
      log_time += time()
      # end objective function

      if (fk <= fk0 - c_1*alpha*SparseAccelerator.dot(dfk, dfk))
        break
      end

      alpha = rho*alpha

      if alpha < 10*eps()
        alpha=0.1
        error("Error in Line search - alpha close to working precision")
      end
    end
    # end backtracking line search using armijo criterion

    #@printf("[%5d]%10.3e%10.3e%10.7f\n", it, alpha, norm(dfk), fk0)

    #x = x + alpha*dx
    SparseAccelerator.WAXPBY!(x, 1, x, alpha, dx)

    spmv_time -= time()
    Xw = SparseAccelerator.SpMV(X, x, fknob_spmv4)
    spmv_count += 1
    spmv_time += time()

    yXw = SparseAccelerator.element_wise_multiply(y, Xw)

    log_time -= time()
    #temp = y./(1+exp(yXw))
    SparseAccelerator.exp!(temp, yXw)
    SparseAccelerator.WAXPB!(temp, 1, temp, 1)
    SparseAccelerator.element_wise_divide!(temp, y, temp)
    log_time += time()

    spmv_time -= time()
    #dfkp1 = -(Xt*temp)/m + lambda*x
    dfkp1 = SparseAccelerator.SpMV(-1/m, Xt, temp, lambda, x, fknob_spmv5)
    spmv_count += 1
    spmv_time += time()

    #S[:,(it - 1)%k + 1] = alpha*dx
    S[:, (it - 1)%k + 1] = SparseAccelerator.WAXPBY(alpha, dx, 0, dx)
    #Y[:,(it - 1)%k + 1] = dfkp1 - dfk
    Y[:, (it - 1)%k + 1] = SparseAccelerator.WAXPBY(1, dfkp1, -1, dfk)
  end

  bw = (nnz(X)*12. + (size(X,1) + size(X,2))*8)*spmv_count/spmv_time/1e9
  println("\nSpMV takes $spmv_time sec ($bw gbps).")
  println("log takes $log_time sec.")
  println("direction takes $direction_time sec.")
  println("lbfgs_opt takes $(time() - t0) sec.")
  x, it
end

function lbfgs_opt_with_reordering(X, y, lambda, xinit, tol, k)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  spmv_time = 0
  log_time = 0
  direction_time = 0
  reorder_time = 0

  m, n = size(X)
  x = xinit
  Xt = X'

  mknobX = (SparseAccelerator.new_matrix_knob)(X, true, true, false, false, false, false)
  mknobXT = (SparseAccelerator.new_matrix_knob)(Xt, true, true, false, false, false, false)

  (SparseAccelerator.set_derivative)(mknobX, SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE, mknobXT)

  fknob_spmv1 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spmv1)
  fknob_spmv2 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobXT, fknob_spmv2)
  fknob_spmv3 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spmv3)
  fknob_spmv4 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spmv4)
  fknob_spmv5 = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobXT, fknob_spmv5)

  (SparseAccelerator.set_reordering_decision_maker)(fknob_spmv1)

  S = zeros(n, k)
  Y = zeros(n, k)

  a = zeros(k, 1)

  dx    = Array(Cdouble, n)

  t0 = time()
  it = 1

  spmv_count = 0

  reordering_status = [false, C_NULL, C_NULL, C_NULL, C_NULL, reorder_time]
  for it=1:100
    spmv_time -= time()
    #Xw = X*x
    Xw = SparseAccelerator.SpMV(X, x, fknob_spmv1)
    spmv_count += 1
    spmv_time += time()

    SparseAccelerator.reordering(
      fknob_spmv1,
      reordering_status,
      Xt, SparseAccelerator.COL_PERM, SparseAccelerator.ROW_INV_PERM,
      :__delimitor__,
      y, SparseAccelerator.ROW_PERM
    )

    #yXw = y.*Xw
    yXw = SparseAccelerator.element_wise_multiply(y, Xw)

    log_time -= time()
    #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
    temp = zeros(length(yXw))
    SparseAccelerator.abs!(temp, yXw)
    SparseAccelerator.WAXPBY!(temp, -1, temp, 0, temp)
    SparseAccelerator.exp!(temp, temp)
    SparseAccelerator.log1p!(temp, temp)
    temp2 = zeros(length(yXw))
    SparseAccelerator.min!(temp2, yXw, 0)
    SparseAccelerator.WAXPBY!(temp, 1, temp, 1, temp2)
    s = SparseAccelerator.sum(temp)

    fk0 = s/m+(lambda/2)*SparseAccelerator.dot(x, x)
    #temp = y./(1+exp(yXw))
    SparseAccelerator.exp!(temp, yXw)
    SparseAccelerator.WAXPB!(temp, 1, temp, 1)
    SparseAccelerator.element_wise_divide!(temp, y, temp)

    log_time += time()

    spmv_time -= time()
    #dfk = -(Xt*temp)/m + lambda*x
    dfk = SparseAccelerator.SpMV(-1/m, Xt, temp, lambda, x, fknob_spmv2)
    spmv_count += 1
    spmv_time += time()

    if (SparseAccelerator.norm(dfk) < tol)
      break;
    end

    direction_time -= time()
    SparseAccelerator.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)
    direction_time += time()

    # backtracking line search using armijo criterion
    alphaMax = 1 # this is the maximum step length
    alpha = alphaMax
    rho = 1/2 # < 1 reduction factor of alpha
    c_1 = 1e-4

    while true
      # logistic loss objective funcion
      #w = x - alpha*dfk
      w = SparseAccelerator.WAXPBY(1, x, -alpha, dfk)

      spmv_time -= time()
      #Xw = X*w
      Xw = SparseAccelerator.SpMV(X, w, fknob_spmv3)
      spmv_count += 1
      spmv_time += time()
      SparseAccelerator.element_wise_multiply!(yXw, y, Xw)

      log_time -= time()
      #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
      SparseAccelerator.abs!(temp, yXw)
      SparseAccelerator.WAXPBY!(temp, -1, temp, 0, temp)
      SparseAccelerator.exp!(temp, temp)
      SparseAccelerator.log1p!(temp, temp)
      temp2 = zeros(length(yXw))
      SparseAccelerator.min!(temp2, yXw, 0)
      SparseAccelerator.WAXPBY!(temp, 1, temp, 1, temp2)
      s = SparseAccelerator.sum(temp)

      fk = s/m+(lambda/2)*SparseAccelerator.dot(w, w)
      log_time += time()
      # end objective function

      if (fk <= fk0 - c_1*alpha*SparseAccelerator.dot(dfk, dfk))
        break
      end

      alpha = rho*alpha

      if alpha < 10*eps()
        alpha=0.1
        error("Error in Line search - alpha close to working precision")
      end
    end
    # end backtracking line search using armijo criterion

    #@printf("[%5d]%10.3e%10.3e%10.7f\n", it, alpha, norm(dfk), fk0)

    #x = x + alpha*dx
    SparseAccelerator.WAXPBY!(x, 1, x, alpha, dx)

    spmv_time -= time()
    Xw = SparseAccelerator.SpMV(X, x, fknob_spmv4)
    spmv_count += 1
    spmv_time += time()

    yXw = SparseAccelerator.element_wise_multiply(y, Xw)

    log_time -= time()
    #temp = y./(1+exp(yXw))
    SparseAccelerator.exp!(temp, yXw)
    SparseAccelerator.WAXPB!(temp, 1, temp, 1)
    SparseAccelerator.element_wise_divide!(temp, y, temp)
    log_time += time()

    spmv_time -= time()
    #dfkp1 = -(Xt*temp)/m + lambda*x
    dfkp1 = SparseAccelerator.SpMV(-1/m, Xt, temp, lambda, x, fknob_spmv5)
    spmv_count += 1
    spmv_time += time()

    #S[:,(it - 1)%k + 1] = alpha*dx
    S[:, (it - 1)%k + 1] = SparseAccelerator.WAXPBY(alpha, dx, 0, dx)
    #Y[:,(it - 1)%k + 1] = dfkp1 - dfk
    Y[:, (it - 1)%k + 1] = SparseAccelerator.WAXPBY(1, dfkp1, -1, dfk)
  end

  SparseAccelerator.reverse_reordering(
    reordering_status,
    :__delimitor__,
    x, SparseAccelerator.COL_PERM)

  bw = (nnz(X)*12. + (size(X,1) + size(X,2))*8)*spmv_count/spmv_time/1e9
  println("\nSpMV takes $spmv_time sec ($bw gbps)")
  println("log takes $log_time sec.")
  println("direction takes $direction_time sec.")
  println("lbfgs_opt takes $(time() - t0) sec.")
  x, it
end


# Load data
println("Loading Data")

filename, fileext = splitext(ARGS[1])
X = matrix_market_read(ARGS[1])
X = SparseMatrixCSC{Float64, Int32}([sparse(ones(size(X, 2), 1)) X'])
y = matrix_market_read(string(filename, "_b", fileext))
y=max(y,0)
(n,p) = size(X)

# Set up problem
lambda = 0
sparsity = nnz(X)/(n*p)
println("n=$n p=$p nnz=$(nnz(X)) stored as sparse=$sparsity")

Xt=X'

w, it = GD(X,Xt,y,lambda, zeros(p), 1e-10)
@printf("Grad-Decent: %d iterations f = %.14f\n", it, LogisticLoss(w,X,Xt,y,lambda)[1])
# Expected output: Grad-Decent: 100 iterations f = 0.34398484995673

w, it = lbfgs_ref(X, y, lambda, zeros(p), 1e-10, 3)
w, it = lbfgs_ref(X, y, lambda, zeros(p), 1e-10, 3)
@printf("Original L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3)
w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3)
@printf("Opt L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

w, it = lbfgs_opt_with_reordering(X, y, lambda, zeros(p), 1e-10, 3)
w, it = lbfgs_opt_with_reordering(X, y, lambda, zeros(p), 1e-10, 3)
@printf("Opt_with_reordering L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
# Expected output: L-BFGS: 33 iterations f = 0.33390367349181

xinit, tol, k = zeros(p), 1e-10, 3
@acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
@printf("Accelerated L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
