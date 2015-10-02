function [x,y,s,f] = pdip(Af) %,bf,pf)
% primal-dual interior-point method for problem % % min p'x s.t. Ax=b, x>=0, % % whose dual is % % max b'y  s.t. A'y+s=p,  s>=0.
%
% calling sequence:
%
% [x,y,s,f] = pdip(A,b,p)
%
% input: A is an m x n SPARSE constraint matrix.
%        b is an m x 1 right-hand side vector
%        p is an n x 1 cost vector.
%
% output: x is the  n x 1 solution of the primal problem
%         y is the m x 1 dual solution
%         s is the n x 1 vector of "dual slacks"
%         f is the optimal objective value

fprintf('%s\n', Af);
% min 2*x1 + x2 subject to x1 + x2 = 1, x1 >= 0, x2 >= 0
% expected solution: x1 = 0, x2 = 1, obj = 1
%A = sparse([1 1]);
%b = [ 1 ]';
%p = [ 2 1 ]';

% small example from http://cvxopt.org/examples/tutorial/lp.html
%A = sparse([ -1 1 1 0 0 ; -1 -1 0 1 0 ; 1 -2 0 0 1 ]);
%b = [ 1 -2 4 ]';
%p = [ 2 1 0 0 0 ]';

A = mmread([Af '-A.mtx']);
b = mmread([Af '-b.mtx']);
p = mmread([Af '-p.mtx']);
if size(p, 1) == 0
  p = ones(size(A,2),1);
end

fprintf('size(A): %d x %d\n', size(A));
if nargin ~= 1
  error('must have three input arguments'); end

if ~issparse(A)
  error('first input argument A must be a SPARSE matrix; possibly use sparse() to convert'); end

t0=cputime;
[m,n] = size(A);
if m <= 0 | n <= 0
  error('input matrix A must be nontrivial'); end

if n ~= length(p)
  error('size of vector p must match number of columns in A'); end
if m ~= length(b)
  error('size of vector b must match number of rows in A'); end

% set initial point, based on largest element in (A,b,p)
bigM = max(max(abs(A))); bigM = max([norm(b,inf), norm(p,inf), bigM]);
x = 100*bigM*ones(n,1); s = x;
y = zeros(m,1);

% find row/column ordering that gives a sparse Cholesky % factorization of ADA'
ordering = symamd(A*A');
norm(b)
norm(p)
size(p)
max(p)
bc = 1+max([norm(b), norm(p)])

for iter=1:200
  
% compute residuals
  Rd = A'*y+s-p;
  Rp = A*x-b;
  Rc = x.*s;
  mu = mean(Rc);
  relResidual = norm([Rd;Rp;Rc])/bc;
%  fprintf('iter %2i: mu = %9.2e, resid = %9.2e\n', iter, mu, relResidual);
  fprintf('iter %2i: log10(mu) = %5.2f, resid = %9.2e, obj = %9.2e\n', iter, log10(full(mu)), ...
          full(relResidual), p'*x);
  if(relResidual <= 1.e-7 & mu <= 1.e-7) break; end;
  Rc = Rc - min(0.1,100*mu)*mu;
  
  % set up the scaling matrix, and form the coefficient matrix for
  % the linear system
  d = min(5.e+15, x./s);
  B = A*sparse(1:n,1:n,d)*A';
  % use the form of the Cholesky routine "cholinc" that's best
  % suited to interior-point methods
  %R = cholinc(B(ordering,ordering),'inf'); % cholin obsolete in recent MATLAB versions
  R = chol(B(ordering, ordering));
  %R = sparse(ldl(full(B(ordering, ordering))))';
  %[L, U] = lu(full(B(ordering, ordering)));
  %R = sparse(diag(sqrt(abs(diag(U))))\U);
  
  % set up the right-hand side
  t1 = x.*Rd-Rc;
  t2 = -(Rp+A*(t1./s));
  
  % solve it and recover the other step components
  dy = zeros(m,1);
  dy(ordering) = R\(R'\t2(ordering));
  dx = (x.*(A'*dy)+t1)./s;
  ds = -(s.*dx+Rc)./x;
  
  tau = max(.9995,1-mu);
  ap = -1/min(min(dx./x),-1);
  ad = -1/min(min(ds./s),-1);
  ap = tau*ap;
  ad = tau*ad;
  x = x + ap*dx;
  s = s + ad*ds;
  y = y + ad*dy;
end

fprintf('iter %2i, resid = %9.2e\n', iter, relResidual)

f = p'*x;

% convert x,y,s to full data structures
x=full(x); s=full(s); y=full(y);

fprintf('Done!\t[m n] = [%g %g]\tCPU = %g\n', m, n, cputime-t0); fprintf('Optimal objective value (f) = %e\n', f); return;  

