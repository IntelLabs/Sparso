#=
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

include("../../src/Sparso.jl")
include("../../src/simple-show.jl")
using Sparso

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS, SA_USE_SPLITTING_PATTERNS)

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

function lbfgs_ref(X, y, lambda, xinit, tol, k)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  Sparso.reset_spmp_spmv_time()
  Sparso.reset_knob_spmv_time()

  bytes = (nnz(X)*12. + (size(X,1) + size(X,2))*8)

  spmv_time = 0
  log_time = 0
  direction_time = 0

  m, n = size(X)
  x = xinit
  Xt = X'

  spmv_count = 0

  S = zeros(n, k)
  Y = zeros(n, k)

  a = zeros(k, 1)

  temp  = Array(Cdouble, m)
  Xw    = Array(Cdouble, m)
  yXw   = Array(Cdouble, m)
  dfk   = Array(Cdouble, n)
  dfkp1 = Array(Cdouble, n) 
  w     = Array(Cdouble, n)
  dx    = Array(Cdouble, n)
  
  t0 = time()
  it = 1
  for it=1:100
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
    Xw = X*x
    spmv_count += 1
    spmv_time += time()

    yXw = y.*Xw

    log_time -= time()
    #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
    #fk0 = s/m+(lambda/2)*norm(x)^2
    fk0 = Sparso.lbfgs_loss_function1(yXw, x, lambda)
    #temp = y./(1+exp(yXw))
    Sparso.lbfgs_loss_function2!(temp, y, yXw)
    log_time += time()

    spmv_time -= time()
    dfk = -(Xt*temp)/m + lambda*x
    spmv_count += 1
    spmv_time += time()

    if (norm(dfk) < tol)
      break;
    end

    # Since Julia code for this is too slow we do this in C to show sufficient speedups by
    # faster SpMV. To be fair, the reference code also uses it.
    direction_time -= time()
    Sparso.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)
    direction_time += time()

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
      spmv_count += 1
      spmv_time += time()
      yXw = y.*Xw

      log_time -= time()
      #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
      #fk = s/m+(lambda/2)*norm(w)^2
      fk = Sparso.lbfgs_loss_function1(yXw, w, lambda)
      log_time += time()
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
    spmv_count += 1
    spmv_time += time()

    yXw = y.*Xw

    log_time -= time()
    Sparso.lbfgs_loss_function2!(temp, y, yXw)
    log_time += time()

    spmv_time -= time()
    dfkp1 = -(Xt*temp)/m + lambda*x
    spmv_count += 1
    spmv_time += time()

    S[:,(it - 1)%k + 1] = alpha*dx
    Y[:,(it - 1)%k + 1] = dfkp1 - dfk
  end

  bytes *= spmv_count
  println("\nSpMV takes $spmv_time sec ($(bytes/spmv_time/1e9) gbps).")

  spmv_time = Sparso.get_spmp_spmv_time()
  if spmv_time > 0
    println("time spent on spmp spmv $spmv_time sec ($(bytes/spmv_time/1e9) gbps)")
    spmv_time = Sparso.get_knob_spmv_time()
    println("time spent on knob spmv $spmv_time sec ($(bytes/spmv_time/1e9) gbps)")
  end

  println("log takes $log_time sec.")
  println("direction takes $direction_time sec.")
  println("lbfgs_ref takes $(time() - t0) sec.")

  x, it
end

function lbfgs_opt(X, y, lambda, xinit, tol, k, do_print, with_context_opt)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  Sparso.reset_spmp_spmv_time()
  Sparso.reset_knob_spmv_time()

  spmv_time = 0
  log_time = 0
  direction_time = 0

  m, n = size(X)
  x = xinit
  Xt = X'

  spmv_count = 0

if with_context_opt
  mknobX = (Sparso.new_matrix_knob)(:X, true, true, false, false, false, false)
  mknobXT = (Sparso.new_matrix_knob)(:Xt, true, true, false, false, false, false)

  (Sparso.set_derivative)(mknobX, Sparso.DERIVATIVE_TYPE_TRANSPOSE, mknobXT)

  fknob_spmv1 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobX, fknob_spmv1)
  fknob_spmv2 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobXT, fknob_spmv2)
  fknob_spmv3 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobX, fknob_spmv3)
  fknob_spmv4 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobX, fknob_spmv4)
  fknob_spmv5 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobXT, fknob_spmv5)
else
  fknob_spmv1 = C_NULL
  fknob_spmv2 = C_NULL
  fknob_spmv3 = C_NULL
  fknob_spmv4 = C_NULL
  fknob_spmv5 = C_NULL
end

  S = zeros(n, k)
  Y = zeros(n, k)

  a = zeros(k, 1)
  dx = zeros(n)
  temp = zeros(m)
  Xw = zeros(m)
  yXw = zeros(m)
  dfk = zeros(n)
  dfk1 = zeros(n)
  dfkp1 = zeros(n)
  w = zeros(n)

  t0 = time()
  it = 1
  for it=1:100
    spmv_time -= time()
    #Xw = X*x
    Sparso.SpMV!(Xw, X, x, fknob_spmv1)
    spmv_count += 1
    spmv_time += time()

    #yXw = y.*Xw
    Sparso.element_wise_multiply!(yXw, y, Xw)

    log_time -= time()
    #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
    #fk0 = s/m+(lambda/2)*norm(x)^2
    fk0 = Sparso.lbfgs_loss_function1(yXw, x, lambda)
    #temp = y./(1+exp(yXw))
    Sparso.lbfgs_loss_function2!(temp, y, yXw)

    log_time += time()

    spmv_time -= time()
    #dfk = -(Xt*temp)/m + lambda*x
    Sparso.SpMV!(dfk, -1/m, Xt, temp, lambda, x, 0, fknob_spmv2)
    spmv_count += 1
    spmv_time += time()

    if (Sparso.norm(dfk) < tol)
      break;
    end

    direction_time -= time()
    Sparso.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)
    direction_time += time()

    # backtracking line search using armijo criterion
    alphaMax = 1 # this is the maximum step length
    alpha = alphaMax
    rho = 1/2 # < 1 reduction factor of alpha
    c_1 = 1e-4

    while true
      # logistic loss objective funcion
      #w = x - alpha*dfk
      Sparso.WAXPBY!(w, 1, x, -alpha, dfk)

      spmv_time -= time()
      #Xw = X*w
      Sparso.SpMV!(Xw, X, w, fknob_spmv3)
      spmv_count += 1
      spmv_time += time()
      Sparso.element_wise_multiply!(yXw, y, Xw)

      log_time -= time()
      #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
      #fk = s/m+(lambda/2)*norm(w)^2
      fk = Sparso.lbfgs_loss_function1(yXw, w, lambda)
      log_time += time()
      # end objective function

      if (fk <= fk0 - c_1*alpha*Sparso.dot(dfk, dfk))
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
    Sparso.WAXPBY!(x, 1, x, alpha, dx)

    spmv_time -= time()
    Sparso.SpMV!(Xw, X, x, fknob_spmv4)
    spmv_count += 1
    spmv_time += time()

    Sparso.element_wise_multiply!(yXw, y, Xw)

    log_time -= time()
    #temp = y./(1+exp(yXw))
    Sparso.lbfgs_loss_function2!(temp, y, yXw)
    log_time += time()

    spmv_time -= time()
    #dfkp1 = -(Xt*temp)/m + lambda*x
    Sparso.SpMV!(dfkp1, -1/m, Xt, temp, lambda, x, 0, fknob_spmv5)
    spmv_count += 1
    spmv_time += time()

    #S[:,(it - 1)%k + 1] = alpha*dx
    S[:, (it - 1)%k + 1] = Sparso.WAXPBY(alpha, dx, 0, dx)
    #Y[:,(it - 1)%k + 1] = dfkp1 - dfk
    Y[:, (it - 1)%k + 1] = Sparso.WAXPBY(1, dfkp1, -1, dfk)
  end

  if do_print
    bytes = (nnz(X)*12. + (size(X,1) + size(X,2))*8)*spmv_count

    println("\nSpMV takes $spmv_time sec ($(bytes/spmv_time/1e9) gbps).")

    spmv_time = Sparso.get_spmp_spmv_time()
    println("time spent on spmp spmv $spmv_time sec ($(bytes/spmv_time/1e9) gbps)")
    spmv_time = Sparso.get_knob_spmv_time()
    println("time spent on knob spmv $spmv_time sec ($(bytes/spmv_time/1e9) gbps)")

    println("log takes $log_time sec.")
    println("direction takes $direction_time sec.")
    println("lbfgs_opt takes $(time() - t0) sec.")
  end
  x, it
end

function lbfgs_opt_with_reordering(X, y, lambda, xinit, tol, k, do_print)
  #@printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")

  Sparso.reset_spmp_spmv_time()
  Sparso.reset_knob_spmv_time()

  spmv_time = 0
  log_time = 0
  direction_time = 0
  reorder_time = 0

  m, n = size(X)
  x = xinit
  Xt = X'

  mknobX = (Sparso.new_matrix_knob)(:X, true, true, false, false, false, false)
  mknobXT = (Sparso.new_matrix_knob)(:Xt, true, true, false, false, false, false)

  (Sparso.set_derivative)(mknobX, Sparso.DERIVATIVE_TYPE_TRANSPOSE, mknobXT)

  fknob_spmv1 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobX, fknob_spmv1)
  fknob_spmv2 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobXT, fknob_spmv2)
  fknob_spmv3 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobX, fknob_spmv3)
  fknob_spmv4 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobX, fknob_spmv4)
  fknob_spmv5 = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(mknobXT, fknob_spmv5)

  (Sparso.set_reordering_decision_maker)(fknob_spmv1)

  spmv_count = 0

  S = zeros(n, k)
  Y = zeros(n, k)

  a = zeros(k, 1)
  dx = zeros(n)
  temp = zeros(m)
  Xw = zeros(m)
  yXw = zeros(m)
  dfk = zeros(n)
  dfk1 = zeros(n)
  dfkp1 = zeros(n)
  w = zeros(n)

  t0 = time()
  it = 1

  reordering_status = [false, C_NULL, C_NULL, C_NULL, C_NULL, reorder_time]
  for it=1:100
    spmv_time -= time()
    #Xw = X*x
    Sparso.SpMV!(Xw, X, x, fknob_spmv1)
    spmv_count += 1
    spmv_time += time()

    Sparso.reordering(
      fknob_spmv1,
      reordering_status,
      Xt, Sparso.COL_PERM, Sparso.ROW_INV_PERM, mknobXT,
      :__delimitor__,
      y, Sparso.ROW_PERM
    )

    #yXw = y.*Xw
    Sparso.element_wise_multiply!(yXw, y, Xw)

    log_time -= time()
    #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
    #fk0 = s/m+(lambda/2)*norm(x)^2
    fk0 = Sparso.lbfgs_loss_function1(yXw, x, lambda)
    #temp = y./(1+exp(yXw))
    Sparso.lbfgs_loss_function2!(temp, y, yXw)
    log_time += time()

    spmv_time -= time()
    #dfk = -(Xt*temp)/m + lambda*x
    Sparso.SpMV!(dfk, -1/m, Xt, temp, lambda, x, 0, fknob_spmv2)
    spmv_count += 1
    spmv_time += time()

    if (Sparso.norm(dfk) < tol)
      break;
    end

    direction_time -= time()
    Sparso.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)
    direction_time += time()

    # backtracking line search using armijo criterion
    alphaMax = 1 # this is the maximum step length
    alpha = alphaMax
    rho = 1/2 # < 1 reduction factor of alpha
    c_1 = 1e-4

    while true
      # logistic loss objective funcion
      #w = x - alpha*dfk
      Sparso.WAXPBY!(w, 1, x, -alpha, dfk)

      spmv_time -= time()
      #Xw = X*w
      Sparso.SpMV!(Xw, X, w, fknob_spmv3)
      spmv_count += 1
      spmv_time += time()
      Sparso.element_wise_multiply!(yXw, y, Xw)

      log_time -= time()
      #s = sum(log(1 + exp(-abs(yXw))) - min(yXw, 0))
      #fk = s/m+(lambda/2)*norm(w)^2
      fk = Sparso.lbfgs_loss_function1(yXw, w, lambda)
      log_time += time()
      # end objective function

      if (fk <= fk0 - c_1*alpha*Sparso.dot(dfk, dfk))
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
    Sparso.WAXPBY!(x, 1, x, alpha, dx)

    spmv_time -= time()
    Sparso.SpMV!(Xw, X, x, fknob_spmv4)
    spmv_count += 1
    spmv_time += time()

    Sparso.element_wise_multiply!(yXw, y, Xw)

    log_time -= time()
    #temp = y./(1+exp(yXw))
    Sparso.lbfgs_loss_function2!(temp, y, yXw)
    log_time += time()

    spmv_time -= time()
    #dfkp1 = -(Xt*temp)/m + lambda*x
    Sparso.SpMV!(dfkp1, -1/m, Xt, temp, lambda, x, 0, fknob_spmv5)
    spmv_count += 1
    spmv_time += time()

    #S[:,(it - 1)%k + 1] = alpha*dx
    S[:, (it - 1)%k + 1] = Sparso.WAXPBY(alpha, dx, 0, dx)
    #Y[:,(it - 1)%k + 1] = dfkp1 - dfk
    Y[:, (it - 1)%k + 1] = Sparso.WAXPBY(1, dfkp1, -1, dfk)
  end

  Sparso.reverse_reordering(
    reordering_status,
    :__delimitor__,
    x, Sparso.COL_PERM)

  if do_print
    bytes = (nnz(X)*12. + (size(X,1) + size(X,2))*8)*spmv_count

    println("\nSpMV takes $spmv_time sec ($(bytes/spmv_time/1e9) gbps).")

    spmv_time = Sparso.get_spmp_spmv_time()
    println("time spent on spmp spmv $spmv_time sec ($(bytes/spmv_time/1e9) gbps)")
    spmv_time = Sparso.get_knob_spmv_time()
    println("time spent on knob spmv $spmv_time sec ($(bytes/spmv_time/1e9) gbps)")

    println("log takes $log_time sec.")
    println("direction takes $direction_time sec.")
    println("lbfgs_opt takes $(time() - t0) sec.")
  end
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

#Xt=X'

#w, it = GD(X,Xt,y,lambda, zeros(p), 1e-10)
#@printf("Grad-Decent: %d iterations f = %.14f\n", it, LogisticLoss(w,X,Xt,y,lambda)[1])
# Expected output: Grad-Decent: 100 iterations f = 0.34398484995673

w, it = lbfgs_ref(X, y, lambda, zeros(p), 1e-10, 3)
w, it = lbfgs_ref(X, y, lambda, zeros(p), 1e-10, 3)
@printf("Original L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, false, false)
w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, true, false)
@printf("Call-repl L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, false, true)
w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, true, true)
@printf("Opt L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

Sparso.set_knob_log_level(1)
w, it = lbfgs_opt_with_reordering(X, y, lambda, zeros(p), 1e-10, 3, false)
w, it = lbfgs_opt_with_reordering(X, y, lambda, zeros(p), 1e-10, 3, true)
Sparso.set_knob_log_level(0)
@printf("Opt_with_reordering L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
# Expected output: L-BFGS: 33 iterations f = 0.33390367349181

xinit, tol, k = zeros(p), 1e-10, 3
@acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
@printf("First accelerated L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
#Sparso.set_knob_log_level(1)
xinit, tol, k = zeros(p), 1e-10, 3
@acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
@printf("Accelerated L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
