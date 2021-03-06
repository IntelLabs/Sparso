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

include("../../../src/Sparso.jl")
using Sparso

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

if length(ARGS) == 2
  test = ARGS[2]
else
  test = "julia"
end

if test == "julia"
  println("compiler warmup (ignored): ")
  w, it = lbfgs_ref(X, y, lambda, zeros(p), 1e-10, 3)

  println("\nRUN: ")
  w, it = lbfgs_ref(X, y, lambda, zeros(p), 1e-10, 3)
else
  if test == "call-repl"
    set_options(SA_ENABLE, SA_USE_SPMP, SA_REPLACE_CALLS)
  elseif test == "context"
    set_options(SA_ENABLE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS)
  elseif test == "reorder"
    set_options(SA_ENABLE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)
  elseif test == "verbose"
    set_options(SA_ENABLE, SA_USE_SPMP, SA_VERBOSE, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)
  end

  println("compiler warnup (ignored): ")
  xinit, tol, k = zeros(p), 1e-10, 3
  @acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)

  println("\nRUN: ")
  xinit, tol, k = zeros(p), 1e-10, 3
  @acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
end

@printf("L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
