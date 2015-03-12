# Pkg.add("MAT") to use MAT package
using MAT

function mylogsumexp(b)
  # does logsumexp across column
  log(1 + exp(-abs(b)))+max(b,0)
end

function GD(X, Xt, y, lambda, xinit)
  x=xinit
  t=0
  @printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")
  for it=1:100
    (n,p) = size(X)
    tic()
    Xw = X*x
    yXw = y.*Xw
    t += toq()

    s = sum(mylogsumexp(-yXw))
    fk0 = s/n+(lambda/2)*norm(x)^2
    dfk = -(Xt*(y./(1+exp(yXw))))/n + lambda*w

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
      yXw = y.*Xw
      t += toq()

      s = sum(mylogsumexp(-yXw))
      fk = s/n+(lambda/2)*norm(w)^2
      # end objective function

      if (fk <= fk0 - c_1*alpha*norm(dfk)^2)
        break
      end

      alpha = rho*alpha

      if alpha < 10*eps()
        alpha=0.1
        return alpha # error('Error in Line search - alpha close to working precision');
      end
    end
    # end backtracking line search using armijo criterion

    @printf("[%5d]%10.3e%10.3e%10.7f\n", it, alpha, norm(dfk), fk0)
    x=x-alpha*dfk
  end
  println("SpMV takes $t sec.")
  println("log takes $t sec.")
  x
end

# Load data
println("Loading Data")
vars = matread("rcv1_train_1padded.binary.mat")
X = vars["X"]
y = vec(vars["y"])
y=max(y,0)
#X = [ones(size(X,1),1) X] # already padded in rcv1_train_1padded. adding 1 column very slow in Julia
(n,p) = size(X)

# Set up problem
lambda = 0
sparsity = nnz(X)/(n*p)
println("n=$n p=$p nnz=$(nnz(X)) stored as sparse=$sparsity")

Xt=X'

@time w = GD(X,Xt,y,lambda, zeros(p))
@printf("Grad-Decen: f = %.14f\n", LogisticLoss(w,X,Xt,y,lambda)[1])
# Expected output: Grad-Decen: f = 0.34398484995673
