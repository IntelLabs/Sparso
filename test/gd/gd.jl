# Pkg.add("MAT") to use MAT package
using MAT

function mylogsumexp(b)
  # does logsumexp across column
  B = maximum(b, 2)
  log(sum(exp(b-repmat(B,1,size(b,2))),2))+B
end

function LogisticLoss2(w, X, Xt, y, lambda)
  (n,p) = size(X)
  Xw = X*w
  yXw = y.*Xw

  fk = 1/n*sum(mylogsumexp([zeros(n,1) -yXw]))+(lambda/2)*norm(w)*norm(w)
  gk = -(Xt*(y./(1+exp(yXw))))/n + lambda*w
  (fk, gk)
end

function GD(lossfunc, xinit, options)
  x=xinit
  @printf("%5s%10s%10s%10s\n", "Iter", "alpha", "|dfk|", "fk")
  for it=1:options["MaxIter"]
    (fk,dfk)=lossfunc(x)
    alpha = backtrackingLineSearch(lossfunc,fk,x,dfk,-dfk)
    @printf("[%5d]%10.3e%10.3e%10.7f\n", it, alpha, norm(dfk), fk)
    x=x-alpha*dfk
  end
  x
end

function backtrackingLineSearch(objFunc,objFuncValue,x,dx,dir)
  # backtracking line search using armijo criterion
  # objFunc      - handle for objective function
  # objFuncValue - current objective function value @ x
  # x            - x
  # dx           - dx
  # dir          - search direction
  #
  # example : mb_backtrackingLineSearch(objFunc,objFuncValue,x,dx,dir)

  alphaMax = 1 # this is the maximum step length
  alpha = alphaMax
  rho = 1/2 # < 1 reduction factor of alpha
  c_1 = 1e-4

  while objFunc(x+alpha*dir)[1][1] > objFuncValue + (c_1*alpha*dir'*dx)[1]
    alpha = rho*alpha

    if alpha < 10*eps
      alpha=0.1
      return alpha # error('Error in Line search - alpha close to working precision');
    end
  end

  alpha
end

# Load data
println("Loading Data")
vars = matread("rcv1_train_1padded.binary.mat")
X = vars["X"]
y = vars["y"]
y=max(y,0)
#X = [ones(size(X,1),1) X] # already padded in rcv1_train_1padded. adding 1 column very slow in Julia
(n,p) = size(X)

# Set up problem
lambda = 0
sparsity = nnz(X)/(n*p)
println("n=$n p=$p nnz=$(nnz(X)) stored as sparse=$sparsity")

Xt=X'
objective(w) = LogisticLoss2(w,X,Xt,y,lambda)

options={"MaxIter" => 100}
w = GD(objective, zeros(p), options)
@printf("Grad-Decen: f = %.14f\n", objective(w)[1])

