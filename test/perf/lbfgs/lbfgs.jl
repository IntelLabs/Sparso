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
    Sparso.lbfgs_compute_direction!(dx, k, it, n, S, Y, dfk)

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

#w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, false)
#w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, false)
#@printf("Call-repl L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

#w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, true)
#w, it = lbfgs_opt(X, y, lambda, zeros(p), 1e-10, 3, true)
#@printf("Opt L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])

#w, it = lbfgs_opt_with_reordering(X, y, lambda, zeros(p), 1e-10, 3)
#w, it = lbfgs_opt_with_reordering(X, y, lambda, zeros(p), 1e-10, 3)
#@printf("Opt_with_reordering L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
# Expected output: L-BFGS: 33 iterations f = 0.33390367349181

original_X = copy(X)
original_y = copy(y)
xinit, tol, k = zeros(p), 1e-10, 3
@acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
@printf("First accelerated L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
# It seems that sparse accelerator generated code has updated xinit internally. That caused
# strange behavior when lbfgs_ref is called again: it would return after only 1 iteration.

X = original_X
y = original_y
xinit, tol, k = zeros(p), 1e-10, 3
@acc w, it = lbfgs_ref(X, y, lambda, xinit, tol, k)
@printf("Accelerated L-BFGS:      %d iterations f = %.14f\n", it, LogisticLoss(w,X,X',y,lambda)[1])
