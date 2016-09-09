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

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

global mainloop_time = 0.
global deflate_time = 0.
global blanz_time = 0.
global eig_time = 0.
global sort_time = 0.
global residual_time = 0.
global convtests_time = 0.
global applyshifts_time = 0.

global blanz_ortho_time = 0.
global blanz_spmv_time = 0.
global blanz_ortho_inner_time = 0.
global blanz_qr_time = 0.
global blanz_etc1_time = 0.
global blanz_qrsingblk_time = 0.
global blanz_etc2_time = 0.

global dgemm1_time = 0.
global dgemm2_time = 0.
global dgemm3_time = 0.
global dgemm4_time = 0.
global dgemm5_time = 0.

const LIB_CSR_PATH = "../../lib/libcsr.so"

function mkl_dgemm!(transA::Bool, transB::Bool, alpha, A::StridedMatrix{Float64}, B::StridedMatrix{Float64}, beta, C::StridedMatrix{Float64})

   const CBLAS_ROW_MAJOR = 101;
   const CBLAS_COL_MAJOR = 102;
   const CBLAS_NO_TRANS = 111;
   const CBLAS_TRANS = 112;

   m = size(A, transA ? 2 : 1)
   k = size(A, transA ? 1 : 2)
   n = size(B, transB ? 1 : 2)

   ccall((:cblas_dgemm, LIB_CSR_PATH), Void,
        (Cint, Cint, Cint,
         Cint, Cint, Cint,
         Cdouble, Ptr{Cdouble}, Cint,
         Ptr{Cdouble}, Cint,
         Cdouble, Ptr{Cdouble}, Cint),
         CBLAS_COL_MAJOR,
         transA ? CBLAS_TRANS : CBLAS_NO_TRANS,
         transB ? CBLAS_TRANS : CBLAS_NO_TRANS,
         m, n, k,
         alpha, A, stride(A, 2),
         B, stride(B, 2),
         beta, C, stride(C, 2))
end

function spalloc(m,n,nzmax)
   SparseMatrixCSC(m, n, Array(Int64, n+1), Array(Int64, nzmax), Array(Float64, nzmax))
end

function orth(X)
  if (size(X,2)>1)
    return(svd(X)[1]);
  else
    return(X/L2(X));
  end
end

function dump_matrix(A, name, file_name)
   fp = open(file_name, "a+")
   @printf(fp, "%s\n", name);
   for i = 1:size(A,1)
      for j = 1:size(A,2)
         @printf(fp, "%.17e ", A[i,j])
      end
      @printf(fp, "\n")
   end
   close(fp)
end

function irbleigs(A)
# IRBLEIGS: Finds a few eigenvalues and eigenvectors of a Hermitian matrix.
#
# IRBLEIGS will find a few eigenvalues and eigenvectors for either the 
# standard eigenvalue problem A*x = lambda*x or the generalized eigenvalue 
# problem A*x = lambda*M*x, where A is a sparse Hermitian matrix and the 
# matrix M is positive definite.
#
# [V,D,PRGINF] = IRBLEIGS(A,OPTIONS) 
# [V,D,PRGINF] = IRBLEIGS(A,M,OPTIONS)
# [V,D,PRGINF] = IRBLEIGS('Afunc',N,OPTIONS) 
# [V,D,PRGINF] = IRBLEIGS('Afunc',N,M,OPTIONS)
#
# The first input argument into IRBLEIGS can be a numeric matrix A or an M-file 
# ('Afunc') that computes the product A*X, where X is a (N x blsz) matrix. If A is 
# passed as an M-file then the (N x blsz) matrix X is the first input argument, the 
# second input argument is N, (the size of the matrix A), and the third input argument 
# is blsz, i.e. Afunc(X,N,blsz). For the generalized eigenvalue problem the matrix M is 
# positive definite and passed only as a numeric matrix. IRBLEIGS will compute the 
# Cholesky factorization of the matrix M by calling the internal MATLAB function CHOL. 
# M may be a dense or sparse matrix. In all the implementations IRBLEIGS(A,...) can
# be replaced with IRBLEIGS('Afunc',N,...).
#
# NOTE: If the M-file 'Afunc' requires additional input parameters, e.g. a data structure, 
#       use the option FUNPAR to pass any additional parameters to your function (X,N,blsz,FUNPAR).
# 
# OUTPUT OPTIONS:
# ---------------
#
# I.)   IRBLEIGS(A) or IRBLEIGS(A,M)
#       Displays the desired eigenvalues.    
#
# II.)  D = IRBLEIGS(A) or D = IRBLEIGS(A,M)
#       Returns the desired eigenvalues in the vector D. 
#
# III.) [V,D] = IRBLEIGS(A) or [V,D] = IRBLEIGS(A,M)  
#       D is a diagonal matrix that contains the desired eigenvalues along the 
#       diagonal and the matrix V contains the corresponding eigenvectors, such 
#       that A*V = V*D or A*V = M*V*D. If IRBLEIGS reaches the maximum number of
#       iterations before convergence then V = [] and D = []. Use output option
#       IV.) to get approximations to the Ritz pairs.
#
# IV.)  [V,D,PRGINF] = IRBLEIGS(A) or [V,D,PRGINF] = IRBLEIGS(A,M)
#       This option returns the same as (III) plus a two dimensional array PRGINF 
#       that reports if the algorithm converges and the number of matrix vector 
#       products. If PRGINF(1) = 0 then this implies normal return: all eigenvalues have 
#       converged. If PRGINF(1) = 1 then the maximum number of iterations have been 
#       reached before all desired eigenvalues have converged. PRGINF(2) contains the 
#       number of matrix vector products used by the code. If the maximum number of 
#       iterations are reached (PRGINF(1) = 1), then the matrices V and D contain any 
#       eigenpairs that have converged plus the last Ritz pair approximation for the 
#       eigenpairs that have not converged.
#
# INPUT OPTIONS:
# --------------
#                                   
#       ... = IRBLEIGS(A,OPTS) or  ... = IRBLEIGS(A,M,OPTS)
#       OPTS is a structure containing input parameters. The input parameters can
#       be given in any order. The structure OPTS may contain some or all of the 
#       following input parameters. The string for the input parameters can contain
#       upper or lower case characters. 
#       
#  INPUT PARAMETER      DESCRIPTION                     
#   
#  OPTS.BLSZ         Block size of the Lanczos tridiagonal matrix.        
#                    DEFAULT VALUE    BLSZ = 3
#                        
#  OPTS.CHOLM        Indicates if the Cholesky factorization of the matrix M is available. If
#                    the Cholesky factorization matrix R is available then set CHOLM = 1 and 
#                    replace the input matrix M with R where M = R'*R.
#                    DEFAULT VALUE   CHOLM = 0
#
#  OPTS.PERMM        Permutation vector for the Cholesky factorization of M(PERMM,PERMM). 
#                    When the input matrix M is replaced with R where M(PERMM,PERMM)=R'*R
#                    then the vector PERMM is the permutation vector. 
#                    DEFAULT VALUE   PERMM=1:N
#
#  OPTS.DISPR        Indicates if K Ritz values and residuals are to be displayed on each 
#                    iteration. Set positive to display the Ritz values and residuals on 
#                    each iteration.
#                    DEFAULT VALUE   DISPR = 0 
#
#  OPTS.EIGVEC       A matrix of converged eigenvectors.        
#                    DEFAULT VALUE  EIGVEC = []
#                                                      
#  OPTS.ENDPT        Three letter string specifying the location of the interior end- 
#                    points for the dampening interval(s). 
#                    'FLT' - Let the interior endpoints float.                  
#                    'MON' - Interior endpoints are chosen so that the size of the 
#                            dampening interval is increasing. This creates a nested
#                            sequence of intervals. The interior endpoint will approach 
#                            the closest Ritz value in the undampened part of the spectrum 
#                            to the dampening interval.
#                    DEFAULT VALUE   ENDPT = 'MON' (If SIGMA = 'LE' or 'SE'.)
#                    DEFAULT VALUE   ENDPT = 'FLT' (If SIGMA = a numeric value NVAL.)
#
#  OPTS.FUNPAR       If A is passed as a M-file then FUNPAR contains any additional parameters 
#                    that the M-file requires in order to compute the matrix vector product. 
#                    FUNPAR can be passed as any type, numeric, character, data structure, etc. The
#                    M-file must take the input parameters in the following order (X,n,blsz,FUNPAR).
#                    DEFAULT VALUE  FUNPAR =[]                   
#
#  OPTS.K            Number of desired eigenvalues.             
#                    DEFAULT VALUE  K = 3
#
#  OPTS.MAXIT        Maximum number of iterations, i.e. maximum number of block Lanczos restarts.                           
#                    DEFAULT VALUE  MAXIT = 100
#
#  OPTS.MAXDPOL      Numeric value indicating the maximum degree of the dampening 
#                    polynomial allowed.  
#                    DEFAULT VALUE   MAXDPOL = 200     (If SIGMA = 'LE' or 'SE'.)
#                    DEFAULT VALUE   MAXDPOL = N       (If SIGMA = a numeric value NVAL.)
#
#  OPTS.NBLS         Number of blocks in the Lanczos tridiagonal matrix. The program may increase
#                    NBLS to ensure certain requirements in [1] are satisfied. A warning message
#                    will be displayed if NBLS increases.                          
#                    DEFAULT VALUE    NBLS = 3
#
#  OPTS.SIGMA        Two letter string or numeric value specifying the location 
#                    of the desired eigenvalues.            
#                    'SE'  Smallest Real eigenvalues.                
#                    'LE'  Largest Real eigenvalues.                 
#                    NVAL  A numeric value. The program searches for the K closest
#                          eigenvalues to the numeric value NVAL. 
#                    DEFAULT VALUE   SIGMA = 'LE'
#                                                         
#  OPTS.SIZINT       Size of the dampening interval. Value of 1 indicates consecutive
#                    Ritz values are used to determine the endpoints of the dampening
#                    interval. Value of 2 indicates endpoints are chosen from Ritz
#                    values that are seprated by a single Ritz value. A value of 3
#                    indicates endpoints are chosen from Ritz values that are seprated 
#                    by two Ritz values. Etc. The minimum value is 1 and the maximum 
#                    value is (NBLS-1)*BLSZ-K. The program may modify SIZINT without
#                    notification to ensure certain requirements in [1] are satisfied. 
#                    DEFAULT VALUE    SIZINT = 1
#                                 
#  OPTS.TOL          Tolerance used for convergence. Convergence is determined when             
#                    || Ax - lambda*x ||_2 <= TOL*||A||_2. ||A||_2 is approximated by
#                    largest absolute Ritz value.  
#                    DEFAULT VALUE    TOL = 1d-6
#                                                              
#  OPTS.V0           A matrix of starting vectors.       
#                    DEFAULT VALUE  V0 = randn
#
#  OPTS.ZERTYP       Two letter string to indicate which type of zeros to apply.              
#                    'WL' - Weighted fast Leja points. The weight functions are used to help
#                           increase convergence.  
#                    'ML' - Mapped fast Leja points. Fast Leja points are computed on [-2,2] 
#                           and mapped to the dampening interval. This option is not available
#                           when sigma is a numeric value NVAL.
#                    DEFAULT VALUE  ZERTYP = 'ML' (If SIGMA = 'LE' or 'SE'.)
#
# 
#  DATE MODIFIED: 04/20/2004
#  VER:  1.0

#  AUTHORS:
#  James Baglama     University of Rhode Island, E-mail: jbaglama@math.uri.edu
#  Daniela Calvetti  Case Western Reserve University,  E-mail: dxc57@po.cwru.edu
#  Lothar Reichel    Kent State University, E-mail: reichel@mcs.kent.edu
#  
# REFERENCES:
#   1.) "IRBL: An Implicitly Restarted Block Lanczos Method for large-scale Hermitian 
#        eigenproblems", J. Baglama, D. Calvetti, and L. Reichel, SIAM J. Sci. Comput., 
#        in press, 2003.
#   2.) "irbleigs: A MATLAB program for computing a few eigenpairs of a large sparse
#        Hermitian matrix", J. Baglama, D. Calvetti, and L. Reichel, Technical Report
#        submitted for publication (2001).
#   3.) "Dealing With Linear Dependence during the Iterations of the Restarted
#        Block Lanczos Methods", J. Baglama, Num. Algs., 25, (2000) pp. 23-36. 
#   4.) "Fast Leja Points", J. Baglama, D. Calvetti, and L. Reichel, ETNA, 
#        Vol. 7 (1998), pp. 124-140.
#   5.) "Computation of a few close eigenvalues of a large matrix with 
#        application to liquid crystal modeling", J. Baglama, D. Calvetti, 
#        L. Reichel, and A. Ruttan, J. of Comp. Phys., 146 (1998), pp. 203-226.
#   6.) "Iterative Methods for the Computation of a Few Eigenvalues of a Large 
#        Symmetric Matrix", J. Baglama, D. Calvetti, and L. Reichel, BIT, 36 
#        (1996), pp. 400-421.


# Values used in the GUI demo IRBLDEMO. Not needed for command line computation.
#global matrixprod dnbls err output waithan;
#string1 = 'plot(0.5,0.5,''Color'',[0.8 0.8 0.8]);';
#string2 = 'set(gca,''XTick'',[],''YTick'',[],''Visible'',''off'');';
#string3 = 'text(0,0.5,err,''FontSize'',10,''Color'',''r'');';
#string4 = 'axis([0 1 0 1]);';
#output  = strcat(string1,string2,string3,string4);

# Too many output arguments requested.
#if (nargout >= 4) error("ERROR: Too many output arguments."); end

#----------------------------#
# BEGIN: PARSE INPUT VALUES. #
#----------------------------#

# No input arguments, return help.
#if length(vararg) == 0
   #help(irbleigs)
   #return
#end

# Get matrix A. Check type (numeric or character) and dimensions.
#if (isstruct(varargin{1})) 
   #err = "A must be a matrix.";
   #if ~isempty(gcbo), close(waithan); eval(output); end, error(err);
#end
#if ischar(A)
   #if nargin == 1, error("Need N (size of matrix A)."); end  
   #n = varargin{2};
   #if ~isnumeric(n) | length(n) > 2 
      #error("Second argument N must be a numeric value."); 
   #end
#else
   (n,n) = size(A);
   #if any(size(A) ~= n)
      #err = "Matrix A is not square.";
      #if ~isempty(gcbo), close(waithan); eval(output); end, error(err);
   #end
   #if ~isnumeric(A)
      #err = "A must be a numeric matrix.";
      #if ~isempty(gcbo), close(waithan); eval(output); end, error(err);
  #end   
  #if nnz(A) == 0 
     #err = "Matrix A contains all zeros.";
     #if ~isempty(gcbo), close(waithan); eval(output); end, error(err); 
  #end
#end

# Set all input options to default values.
M = []; blsz = 3; cholM = 0; dispr = 0; eigvec = []; endpt = "MON"; K = 3; maxdpol = 200; maxit = 100;  
nbls::Int = 3
nval = []; zertyp = "ML"; sigma = "LE"; sizint = 1; tol = 1e-6; permM = []; funpar = [];

# Set indicators for ENDPTS and MAXDPOL to be empty arrays. The indicators determine which
# default values should be used.
IENDPT = []; IMAXDPOL = [];

# If a generalized eigenvalue problem is to be solved get matrix M. 
# Check type (numeric or character) and dimensions.
#if ((ischar(A) & nargin > 2) | (isnumeric(A) & nargin > 1)) 
   #if ~isstruct(varargin{2+ischar(A)})
      #M = varargin{2+ischar(A)}
      #if ~isnumeric(M)
         #err = "M must be a numeric matrix.";
         #if ~isempty(gcbo)
           #close(waithan); eval(output);
         #end
         #error(err)
      #end   
      #if any(size(M) ~= n)
         #err = "Matrix M must be the same size as A.";
         #if ~isempty(gcbo)
           #close(waithan); eval(output);
         #end
         #error(err)
      #end
   #end
#end   

# Preallocate memory for large matrices.
#V = spalloc(n,nbls*blsz,n*blsz*nbls)
F = spalloc(n,blsz,n*blsz)

# Get input options from the data structure.
K = 5
#if nargin > 1 + ischar(A) + ~isempty(M)
   #options = varargin{2+ischar(A) + ~isempty(M):nargin};
   #names = fieldnames(options);
   #I = strmatch("BLSZ",upper(names),"exact");
   #if ~isempty(I), blsz = getfield(options,names{I}); end
   #I = strmatch("CHOLM",upper(names),"exact");
   #if ~isempty(I), cholM = getfield(options,names{I}); end
   #I = strmatch("DISPR",upper(names),"exact");
   #if ~isempty(I), dispr = getfield(options,names{I}); end
   #I = strmatch("EIGVEC",upper(names),"exact");
   #if ~isempty(I), eigvec = getfield(options,names{I}); end
   #I = strmatch("ENDPT",upper(names),"exact"); IENDPT = I;
   #if ~isempty(I), endpt = upper(getfield(options,names{I})); end
   #I = strmatch("FUNPAR",upper(names),"exact");
   #if  ~isempty(I), funpar = getfield(options,names{I}); end
   #I = strmatch("K",upper(names),"exact");
   #if  ~isempty(I), K = getfield(options,names{I}); end
   #I = strmatch("MAXDPOL",upper(names),"exact"); IMAXDPOL = I;
   #if ~isempty(I), maxdpol = getfield(options,names{I}); end
   #I = strmatch("MAXIT",upper(names),"exact");
   #if ~isempty(I), maxit = getfield(options,names{I}); end
   #I = strmatch("PERMM",upper(names),"exact");
   #if ~isempty(I), permM = getfield(options,names{I}); end
   #I = strmatch("NBLS",upper(names),"exact");
   #if ~isempty(I), nbls = getfield(options,names{I}); end
   #I = strmatch("ZERTYP",upper(names),"exact");
   #if ~isempty(I), zertyp = upper(getfield(options,names{I})); end
   #I = strmatch("SIGMA",upper(names),"exact");
   #if  ~isempty(I), sigma = upper(getfield(options,names{I})); end
   #I = strmatch("SIZINT",upper(names),"exact");
   #if ~isempty(I), sizint = getfield(options,names{I}); end
   #I = strmatch("TOL",upper(names),"exact");
   #if ~isempty(I), tol = getfield(options,names{I}); end
   #I = strmatch("V0",upper(names),"exact");
   #if ~isempty(I), V = getfield(options,names{I}); end
#end 

# If starting matrix V0 is not given then set starting matrix V0 to be a 
# (n x blsz) matrix of normally distributed random numbers.
#if nnz(V) == 0
  V = randn(n,blsz);
#end 

# Check type of input values and output error message if needed.
#if (!isnumeric(blsz) || !isnumeric(cholM)   || !isnumeric(dispr) || !ischar(endpt)   || 
    #!isnumeric(K)    || !isnumeric(maxdpol) || !isnumeric(maxit) || !isnumeric(nbls) || 
    #!ischar(zertyp)   || !isnumeric(sizint)  || !isnumeric(tol) || !isnumeric(permM))
   #error("Incorrect type for input value(s) in the structure.");
#end

# If a numeric value is given for sigma then set nval=sigma and sigma = "IE" to
# denote that the code is searching for interior eigenvalues.
#if isnumeric(sigma)
  #nval = sigma; sigma = "IE";
#end

# Check the length of the character values sigma, endpt, and zertyp.
if length(sigma) != 2
   err = "SIGMA must be SE, LE, or a numeric value"; 
   if !isempty(gcbo)
     close(waithan); eval(output)
   end
   error(err);
end
if length(endpt) != 3
   error("Incorrect value for ENDPT")
end
if length(zertyp) != 2
   error("Incorrect value for ZERTYP")
end
 
# Resize Krylov subspace if blsz*nbls (i.e. number of Lanczos vectors) is larger 
# than n (i.e. the size of the matrix A).
if blsz*nbls >= n
  nbls = Int(floor(n/blsz))
end

# Check for input errors in the data structure.
#if sizint < 1,  error("Incorrect value for SIZINT. SIZINT must be >= 1."), end
#if K     <= 0,  error("K must be a positive value."), end
#if K     >  n,  error("K must be less than the size of A."), end   
#if blsz  <= 0,  error("BLSZ must be a positive value."), end
#if nbls  <= 1,  error("NBLS must be greater than 1."),   end
#if tol   <  0,  error("TOL must be non-negative."),      end
#if maxit <= 0,  error("MAXIT must be positive."),        end
#if ((cholM ~= 0) & (cholM ~= 1)), error("Unknown value for cholM."), end
#if ~isempty(permM)
   #if ((size(permM,1) ~= n) | (size(permM,2) ~= 1)) & ...
      #((size(permM,1) ~= 1) | (size(permM,2) ~= n)), error("Incorrect size for PERMM"), end
#end
#if maxdpol < nbls, error("MAXDPOL must be >= NBLS"),   end
#if blsz*nbls - sizint <= blsz, error("SIZINT is too large"), end
#if blsz*nbls - K - blsz - sizint < 0
   #nbls = ceil((K+sizint+blsz)/blsz+0.1); 
   #warning(["Increasing NBLS to ",num2str(nbls)]);
#end
#if blsz*nbls >= n, error("K or SIZINT are too large."), end
#if (~strcmp(sigma,"SE") & ~strcmp(sigma,"LE") & ~strcmp(sigma,"IE"))
   #err = "SIGMA must be SE, LE or a numeric value.";
   #if ~isempty(gcbo),close(waithan);eval(output); end, error(err);   
#end
#if (~strcmp(endpt,"FLT") & ~strcmp(endpt,"MON")), error("Unknown value for ENDPT."),  end
#if (~strcmp(zertyp,"WL") & ~strcmp(zertyp,"ML")), error("Unknown value for ZERTYP."), end
#if strcmp(sigma,"IE")
   #zertyp = "WL"; 
   #if isempty(IENDPT),   endpt = "FLT"; end
   #if isempty(IMAXDPOL), maxdpol = n;   end
#end
#if ~isempty(eigvec)
   #if ~isnumeric(eigvec), error("Incorrect type for input eigenvector(s)."), end
   #if ((size(eigvec,1) ~= n) | (size(eigvec,2) >= n))
      #error("Incorrect size of eigenvector matrix EIGVEC.");
   #end
#end
#if ~isnumeric(V), error("Incorrect starting matrix V0."), end
#if ((size(V,1) ~= n) | (size(V,2) ~= blsz)), error("Incorrect size of starting matrix V0."), end
 
# Set tolerance to machine precision if tol < eps.
if tol < eps()
  tol = eps()
end

# If a generalized eigenvalue problem is to be solved then get the Cholesky factorization 
# of the matrix M (i.e. M=R"*R) and set M=R the upper triangular matrix.
if !isempty(M)
    
   # If M is sparse then compute the sparse cholesky factorization. Using symmmd or 
   # symamd to get a permutation vector permM such that M(permM,permM) will have a 
   # sparser cholesky factorization than M. The permutation is given by
   # M(permM,permM) = I(:,permM)"*M*I(:,permM), where I is the identity matrix. 
   # To compute the product A(permM,permM)*x use y(permM,:) = x; y=A*y; y=y(permM,:);     
   # Use the default value of permM, if permM is not part of the input structure.   
   if isempty(permM)
     permM = 1:n
   end
 
   # If cholM = 1, then the input matrix M is given as R in the Cholesky factorization.
   if !cholM
      (M,cholerr) = chol(M[permM,permM]); 
      if cholerr
         err = "Matrix M must be positive definite";
         if !isempty(gcbo)
           close(waithan); eval(output)
         end
         error(err);
      end
   end
end   

#--------------------------%
# END: PARSE INPUT VALUES. %
#--------------------------%

#-----------------------------------------------------------%
# BEGIN: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
#-----------------------------------------------------------%

# Initialization and description of local variables.
conv = 0;            # Number of desired eigenvalues that have converged.
deflate = false;         # Value used to determine if an undesired eigenvector has converged.
eigresdisp = [];     # Holds the residual values of converged eigenvalues.
eigval = Float64[];         # Array used to hold the converged desired eigenvalues.
iter = 1;            # Main loop iteration count.
fLejapts = Float64[];       # Stores fast Leja points.
lcandpts = [];       # Stores "left" candidate points.
lindex = [];         # Index for lcandpts.
lprd = [];           # Stores the Leja polynomial product for "left" candidate points.
leftendpt = [];      # Stores the left most endpoint of the dampening interval.
leftintendpt = [];   # Stores the interior left endpoint of the dampening interval.
leftLejapt = [];     # Place holder in fast leja routines.
mprod = 0;           # The number of matrix vector products.
norlpol = [];        # Normalizes the Leja polynomial to avoid underflow/overflow.
numbls = nbls;       # Initial size of the Lanczos blocks.
pritz = [];          # Holds previous iteration of Ritz values. Used to determine stagnation.
rcandpts = Float64[];       # Stores "right" candidate points.
rindex = [];         # Index for rcandpts.
ritzconv = "F";      # Boolean to determine if all desired eigenvalues have converged.
rprd = Float64[];           # Stores the Leja polynomial product for "right" candidate points.
rightendpt = [];     # Stores the right most endpoint of the dampening interval.
rightintendpt = [];  # Stores the interior right endpoint of the dampening interval.
rightLejapt = [];    # Place holder in fast Leja routines.
singblk = [];        # Integer values used to indicate singular block(s) in blanz.
sqrteps = sqrt(eps()); # Square root of machine tolerance used in convergence testing.
flcount = 0;         # Count used to determine when the maximum number of shifts is reached.

# Value use for demo only, holds the change in nbls. Not needed for command line computation.
dnbls = nbls;

# Holds the maximum absolute value of all computed Ritz values.
global Rmax; Rmax = [];  

# Determine if eigenvectors are requested.
#if (nargout > 1) computvec = "T"
#else
   computvec = "F"
#end 

# Determine the number of input eigenvector(s).
if !isempty(eigvec)
  ninpeigvec = size(eigvec,2)
else ninpeigvec = 0
end  

#--------------------------------------------------------------------%
# END: DESCRIPTION AND INITIALIZATION OF LOCAL AND GLOBAL VARIABLES. %
#--------------------------------------------------------------------%

#----------------------------%
# BEGIN: MAIN ITERATION LOOP %
#----------------------------%

global mainloop_time -= time()

while (iter <= maxit)
   
   # Check and deflate the number of Lanczos blocks if possible. Do not
   # deflate if searching for interior eigenvalues or if the users has
   # inputed eigenvectors. This causes nbls to increase.
   global deflate_time -= time()
   if (nbls > 2 && sigma != "IE" && ninpeigvec == 0) 
      nbls = numbls - Int(floor(size(eigvec,2)/blsz));
      if (sizint >= nbls*blsz-(abs(K-size(eigvec,2))))
         nbls = Int(max(floor((sizint + abs(K-size(eigvec,2)))/blsz) + 1,numbls)); 
      end
   end
   deflate_time += time()
      
   # Compute the block Lanczos decomposition.
   global blanz_time -= time()
@acc   (F,T,V,blsz,mprod,singblk) = blanz(A,K,M,permM,V,blsz,eigvec,funpar,mprod,n,nbls,singblk,sqrteps,tol);
   
   # Determine number of blocks and size of the block tridiagonal matrix T.
   Tsz = size(T,1); nbls = Int(Tsz/blsz); 
   if floor(nbls) != nbls
      # Reset the starting matrix V and re-compute the block Lanczos decomposition if needed.
      (F,T,V,blsz,mprod,singblk) = blanz(A,K,M,permM,randn(n,blsz),blsz,eigvec,mprod,n,nbls,singblk,sqrteps,tol);
      Tsz = size(T,1); nbls = Int(Tsz/blsz); 
      if floor(nbls) != nbls 
         err="blanz in irbleigs.m returns an incorrect matrix T."; 
         if !isempty(gcbo)
            close(waithan); eval(output)
         end
         error(err);   
      end 
   end   
   blanz_time += time()
      
   # Compute eigenvalues and eigenvectors of the block tridiagonal matrix T. 
   global eig_time -= time()
   (ritz,ritzvec) = eig(T);
   eig_time += time()
      
   # Sort the eigenvalues from smallest to largest.
   global sort_time -= time()
   J = sortperm(real(ritz))
   ritz = ritz[J]
   ritzvec = ritzvec[:,J]
   sort_time += time()
   
   # Reached maximum number of iterations, exit main loop.
   if iter >= maxit
      break
   end;
     
   # Compute the residuals for all ritz values. 
   global residual_time -= time()
   #Y = F*ritzvec[Tsz-(blsz-1):Tsz,:]
   Y = zeros(size(F,1), size(ritzvec,2))
   mkl_dgemm!(
      false, false,
      1, F, ritzvec[Tsz-(blsz-1):Tsz,:], 0, Y)
   residuals = sqrt(sum(Y.*Y, 1))
   residual_time += time()
    
   # Convergence tests.
   global convtests_time -= time()
   conv,deflate,eigval,eigvec,eigresdisp,pritz,ritzconv,singblk =
   convtests(computvec,conv,deflate,dispr,eigval,eigvec,eigresdisp,
             iter,K,nval,pritz,residuals,ritz,ritzconv,ritzvec,sigma,
             singblk,sqrteps,tol,Tsz,V);
   convtests_time += time()
   
   # If all desired Ritz values converged then exit main loop.
   if ritzconv == "T"
      break
   end
   
   # Determine dampening intervals, Leja points, and apply Leja zeros as shifts.
   global applyshifts_time -= time()
   (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,
   rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex,V,flcount) =
   applyshifts(blsz,endpt,F,fLejapts,iter,K,lcandpts,leftendpt,leftintendpt,
               lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,rcandpts,
               rightendpt,rightintendpt,ritz,ritzvec,rprd,rightLejapt,rindex,
               zertyp,sigma,singblk,sizint,sqrteps,Tsz,T,V,flcount);
   applyshifts_time += time()
      
   # Check to see if an undesired Ritz vector has been chosen to be deflated. If so,
   # reset the endpoints of the dampening interval(s) on the next interation.
   if deflate 
      deflate = false; Rmax = []; leftintendpt = []; rightintendpt = []; leftendpt = []; rightendpt = [];
   end
  
   # Update the main iteration loop count.
   iter = iter + 1;
   
end # While loop. 

mainloop_time += time()

#--------------------------%
# END: MAIN ITERATION LOOP %
#--------------------------%

#-----------------------%
# BEGIN: OUTPUT RESULTS %
#-----------------------%

# Test to see if maximum number of iterations is reached.
# If output options 1 or 2 are used then set eigval and
# eigvec to be empty arrays.
# If output option 4 is used then use the last Ritz
# pairs as estimate for the unconverged eigenpairs.
PRGINF = Int[]
push!(PRGINF, 0)
if iter >= maxit && nargout <= 2
   push!(PRGINF, 1); eigval = []; eigvec = [];
end        
if iter >= maxit && nargout == 3
   push!(PRGINF, 1)
   NC = K - length(eigval);
   if sigma == "SE"
      JI = 1:NC
   end
   if sigma == "LE"
      JI = Tsz-NC+1:Tsz
   end
   if sigma == "IE"
     JI = sortperm(abs(ritz-nval))
   end
   eigval = [eigval;ritz[JI[1:NC]]]
   eigvec = [eigvec V*ritzvec[:,JI[1:NC]]]
end   
  
# Sort output arguments.
if !isempty(eigval)
   K = min(length(eigval),K); eigval = real(eigval);
   if sigma == "SE"
      I = sortperm(eigval)
      I = I[1:K]
   end   
   if sigma == "LE"
      I = sortperm(eigval); I = flipdim(I, 1); I = I[1:K];
   end   
   if sigma == "IE"
      I = sortperm(abs(eigval-nval)); I = I[1:K];
      J = sortperm(eigval[I]); I = I[J];
   end
   eigval = eigval[I];
end  

# Output option I: Display eigenvalues only.
#if (nargout == 0)
   eigenvalues = eigval
#end

# Output option II: Set eigenvalues equal to output vector.  
#if (nargout == 1)
   #varargout = (eigval)
#end 

# Output option III and IV: Output diagonal matrix of eigenvalues and
# corresponding matrix of eigenvectors.
#if ((nargout == 2) || (nargout == 3))
   #if !isempty(eigvec)
      #eigvec = eigvec[:,ninpeigvec+1:size(eigvec,2)];
      #eigvec = eigvec[:,I];
      ## Must solve a linear system to extract generalized eigenvectors.
      #if !isempty(M) 
         #eigvec = M\eigvec[permM,:];
      #end
   #else
      #eigvec = []; 
   #end 
   #varargout = (eigvec, diagm(eigval))
#end

# Output option IV: Output PRGINF. 
#if nargout == 3
   #PRGINF(2) = mprod; varargout{3} = PRGINF;
#end

# Used in the GUI demo IRBLDEMO. Not used for command line computation.
# Outputs number of matrix-vector products.
#matrixprod = mprod
#varargout
return eigenvalues, iter
end

#---------------------#
# END: OUTPUT RESULTS #
#---------------------#

#------------------------------------#
# BEGIN: BLOCK LANCZOS DECOMPOSITION #
#------------------------------------#

function blanz(A,K,R,permM,V,blsz,eigvec,funpar,mprod,n,nbls,singblk,sqrteps,tol) 
# Computes the Block Lanczos decomposition, A*V = V*T + F*E^T
# with full reorthogonalization. If the generalized eigenvalue 
# problem A*x = lambda*M*x is to be solved then M=R"*R and the 
# Lanczos decomposition, inv(R")*A*inv(R)*V = V*T + F*E^T is returned.
#
# The matrix A can be passed as a numeric matrix or as a filename.
# However, the matrix M must be passed as a numeric matrix. Note that
# if the matrix A is a filename then the file must accept 
# [X(n,blsz),n,blsz] or [X(n,blsz),n,blsz,funpar] as input in that 
# order and the file must return the matrix product A*X(n,blsz).

# James Baglama
# DATE: 11/06/01

#dump_matrix(V, "V", "julia_V.log")
# Values used in the GUI demo IRBLDEMO. Not needed for command line computation.
global err, output, waithan

# Initialization of residual matrix F, main loop count J, and integer singblk.
F=zeros(n,blsz); J = 1; singblk = Array(Bool, 1); singblk[1] = false; 

# Check size of input vectors V and eigvec for errors.
if size(V,1) != n || blsz < 1
   err = "Incorrect size of starting vectors, V in blanz.";
   error(err); 
end
if size(V,2) != blsz
   blsz = size(V,2)
end
if !isempty(eigvec) 
   if size(eigvec,1) != n 
      err = "Incorrect size of EIGVEC in blanz."; 
      error(err);
   end
end

global blanz_ortho_time -= time()

# First orthogonalization of starting vectors.
V = orth(V);

# Orthogonalize V against all converged eigenvectors. 
if !isempty(eigvec)
   if size(eigvec,2) < size(V,2)
      #V = V - eigvec*(eigvec'*V); doteV = eigvec'*V;
      doteV = zeros(size(eigvec, 2), size(V, 2))
      mkl_dgemm!(
         true, false,
         1, eigvec, V, 0, doteV)
      mkl_dgemm!(
         false, false,
         -1, eigvec, doteV, 1, V)
      mkl_dgemm!(
         true, false,
         1, eigvec, V, 0, doteV)
   else
      V = V - eigvec*(V'*eigvec)'; doteV = (V'*eigvec)';
   end   
   if norm(doteV) > sqrteps
      V = V - eigvec*doteV
   end 
end

# First check of linear dependence of starting vector(s). If starting vector(s)
# are linearly dependent then add normalized random vectors and reorthogonalize.
if rank(V,sqrteps*n*blsz) < blsz
   V = V + randn(n,blsz); V = orth(V);
   if !isempty(eigvec)
      if size(eigvec,2) < size(V,2)
         #V = V - eigvec*(eigvec'*V); doteV = eigvec'*V;
         doteV = zeros(size(eigvec, 2), size(V, 2))
         mkl_dgemm!(
            true, false,
            1, eigvec, V, 0, doteV)
         mkl_dgemm!(
            false, false,
            -1, eigvec, doteV, 1, V)
         mkl_dgemm!(
            true, false,
            1, eigvec, V, 0, doteV)
      else
         V = V - eigvec*(V'*eigvec)'; doteV = (V'*eigvec)';
      end   
      if norm(doteV) > sqrteps
         V = V - eigvec*doteV
      end
   end
end   

# If needed second orthogonalization of starting vectors.
if !isempty(eigvec)
   V = orth(V)
end

# Second check of linear dependence of starting vector(s). If starting vector(s)
# are linearly dependent then add normalized random vectors and reorthogonalize.
if (size(V,2) < blsz)
   V = [V randn(n,blsz-size(V,2))];   
   if !isempty(eigvec)
      if size(eigvec,2) < size(V,2)
         V = V - eigvec*(eigvec'*V); doteV = eigvec'*V;
      else
         V = V - eigvec*(V'*eigvec)'; doteV = (V'*eigvec)';
      end   
      if norm(doteV) > sqrteps
         V = V - eigvec*doteV
      end
   end
   # If needed third orthogonalization of starting vectors.
   V = orth(V);
   # Third check of linear dependence of starting vector(s). If starting vector(s)
   # are linearly dependent fatal error return.
   if (size(V,2) < blsz)
      err = "Dependent starting vectors in block Lanczos."; 
      error(err); 
   end   
end   

blanz_ortho_time += time()

# Check desired size of Lanczos matrix T. If size is greater than n
# reduce size to the next multiple of blsz that is less than n.
if blsz*nbls >= n
   nbls = Int(floor(n/blsz))
end

# Pre-allocate matrices
T = zeros(nbls*blsz, nbls*blsz)
D = zeros(blsz, blsz)
V = [V zeros(size(V, 1), (nbls - 1)*blsz)]
dotFV = zeros(size(V, 2), blsz)

# Begin of main iteration loop for the block Lanczos decomposition.
while (J <= nbls)
   
   # Values used for indices.
   Jblsz = J*blsz; Jm1blszp1 = blsz*(J-1)+1;

   global blanz_spmv_time -= time()
   
   # Matrix product with vector(s).
   F = A*V[:,Jm1blszp1:Jblsz];

   blanz_spmv_time += time()

   global blanz_ortho_inner_time -= time()
      
   # Count the number of matrix vector products.
   mprod = mprod + blsz
      
   global dgemm1_time -= time()
   # Orthogonalize F against the previous Lanczos vectors.
   if (J > 1)
      #F = F - V[:,blsz*(J-2)+1:Jm1blszp1-1]*T[Jm1blszp1:Jblsz,blsz*(J-2)+1:Jm1blszp1-1]'
      mkl_dgemm!(
         false, true,
         -1, V[:,blsz*(J-2)+1:Jm1blszp1-1],
         T[Jm1blszp1:Jblsz,blsz*(J-2)+1:Jm1blszp1-1],
         1, F)
   end
   dgemm1_time += time()
     
   # Compute the diagonal block of T.
   global dgemm2_time -= time()

   # fat*tall multiplication
   #D = F'*V[:,Jm1blszp1:Jblsz];
   mkl_dgemm!(
      true, false,
      1, F, V[:,Jm1blszp1:Jblsz], 0, D)
   dgemm2_time += time()
          
   # One step of the block classical Gram-Schmidt process. 
   global dgemm3_time -= time()

   #F = F - V[:,Jm1blszp1:Jblsz]*D;
   mkl_dgemm!(
      false, false,
      -1, V[:,Jm1blszp1:Jblsz], D, 1, F)

   dgemm3_time += time()
   
   # Full reorthogonalization step.
   global dgemm4_time -= time()

   #dotFV = (F'*V[:,1:Jblsz])'
   mkl_dgemm!(
      true, false,
      1, V[:,1:Jblsz], F, 0, dotFV[1:Jblsz,:])

   #F = F - V[:,1:Jblsz]*dotFV[1:Jblsz,:]
   mkl_dgemm!(
      false, false,
      -1, V[:,1:Jblsz], dotFV[1:Jblsz,:], 1, F)

   if norm(dotFV)>sqrteps
      println("1")
      dotFV2 = (F'*V[:,1:Jblsz])'; dotFV=dotFV+dotFV2; F = F - V[:,1:Jblsz]*dotFV2
   end 
   for i = 1:J
      D = D + dotFV[blsz*(i-1)+1:blsz*i,:]
   end 
   dgemm4_time += time()
            
   # Orthogonalize F against all converged eigenvectors.
   global dgemm5_time -= time()
   if !isempty(eigvec)
      if size(eigvec,2) < size(F,2)
         #F = F - eigvec*(eigvec'*F);
         doteF = zeros(size(eigvec, 2), size(F, 2))
         mkl_dgemm!(
            true, false,
            1, eigvec, F, 0, doteF)
         mkl_dgemm!(
            false, false,
            -1, eigvec, doteF, 1, F)

         #doteF = eigvec'*F; 
         mkl_dgemm!(
            true, false,
            1, eigvec, F, 0, doteF)
      else
         println("3")
         F = F - eigvec*(F'*eigvec)';
         doteF = (F'*eigvec)';
      end 
      if norm(doteF) > sqrteps
         println("4")
         F = F - eigvec*doteF
      end
   end 
   dgemm5_time += time()
    
   # To ensure a symmetric matrix T is computed.
   #if size(T,1) == Jm1blszp1 - 1 && size(T,2) == Jm1blszp1 - 1
      #T = [T zeros(Jm1blszp1-1, blsz); zeros(blsz, Jm1blszp1-1) (tril(D, -1) + tril(D)')]
   #else
      #T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz] = tril(D,-1) + tril(D)';
   #end
   T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz] = tril(D,-1) + tril(D)';

   blanz_ortho_inner_time += time()
              
   # Compute QR factorization and off diagonal block of T.  
   if (J < nbls) 

      global blanz_qr_time -= time()
      (V[:,Jblsz+1:blsz*(J+1)], T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz]) = qr(F)
      blanz_qr_time += time()

      global blanz_etc1_time -= time()

      #V = [V V_new]
      #T = [T; zeros(blsz, Jm1blszp1-1) T_new]
      #T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz] = T_new
      #(V[:,Jblsz+1:blsz*(J+1)],T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz]) = qr(F);
           
      # Check for linearly dependent vectors among the blocks.
      stol = max(min(sqrteps,tol),eps()*maximum(abs(diag((T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz])))));
      I = find(abs(diag(T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz])) .<= stol); 

      blanz_etc1_time += time()
      
      # Linearly dependent vectors detected in the current block.
      if !isempty(I)
                                               
         # Exit. Convergence or not enough vectors to continue to build up the space.
         if !isempty(eigvec)
            sizevec = size(eigvec,2)
         else sizevec = 0
         end
         if ((size(I,1) == blsz) && (size(T,2) >= K)) || (sizevec + size(T,2) >= n)
         
            # Resize T and V and exit.
            T = T[1:size(T,2),1:size(T,2)]; V = V[:,1:size(T,2)]; return; 
         end
         
         # Full Reorthogonalization step to ensure orthogonal vectors.
         V = V[:,1:Jblsz];
         dotFV = (F'*V)'; F = F - V*dotFV;
         if norm(dotFV) > sqrteps
            dotFV2 = (F'*V)'; dotFV = dotFV+dotFV2; F = F - V*dotFV2
         end 
         for i = 1:J 
            T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz] = T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz] + dotFV[blsz*(i-1)+1:blsz*i,:]
         end
         
         # Orthogonalize F against all converged eigenvectors.
         if !isempty(eigvec)
            if size(eigvec,2) < size(F,2)
               F = F - eigvec*(eigvec'*F);  doteF = eigvec'*F;
            else
               F = F - eigvec*(F'*eigvec)'; doteF = (F'*eigvec)';
            end   
            if norm(doteF) > sqrteps
               F = F - eigvec*doteF
            end
         end     
         
         # To ensure a symmetric matrix T is computed.
         T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz] = tril(T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz],-1)+tril(T[Jm1blszp1:Jblsz,Jm1blszp1:Jblsz])';

         # Re-compute QR with random vectors.
         global blanz_qrsingblk_time -= time()
         (V[:,Jblsz+1:blsz*(J+1)],T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz]) =
            qrsingblk(F,V[:,1:Jblsz],eigvec,I,blsz,sqrteps)
         blanz_qrsingblk_time += time()
         
         # Set the singular block indicator to true along with which vector(s)
         # are linearly dependent.
         singblk[1] = true; singblk = [singblk Jblsz+I']; 
         
      end
      
      # Set off diagonal blocks to be equal.
      global blanz_etc2_time -= time()
      #T = [T [zeros(Jm1blszp1-1,blsz); T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz]'; zeros(blsz, blsz)]]
      T[Jm1blszp1:Jblsz,Jblsz+1:blsz*(J+1)] = T[Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz]';
      blanz_etc2_time += time()
      
   end
       
   # Update iteration count (block Lanczos while loop).         
   J = J + 1; 
    
end # While loop.

   return (F,T,V,blsz,mprod,singblk)
end

function qrsingblk(F,V,eigvec,I,blsz,sqrteps)
# This function computes the  QR factorization for a singular 
# off diagonal block of the block tridiagonal matrix T.
# The diagonal element R(k,k) of R associated with the linearly
# dependent vector F(:,k) is set to zero and F(:,k) is set to a
# random vector that is orthogonal to previous Lanczos vectors
# and any converged eigenvectors.

# James Baglama 
# DATE: 11/06/01

# Initialize off-diagonal block to zero.
R = zeros(blsz,blsz); n = size(F,1);

for k = 1 : blsz
   
   if !all(I-k .!= 0)
      # Set F(:,k) to a random vector and orthogonalizes F[:,k] to
      # all previous Lanczos vectors and any converged eigenvectors.
      F[:,k] = randn(n,1); 
      dotVF = V'*F[:,k]; F[:,k] = F[:,k] - V*dotVF; 
      if k > 1
         dotVF = F[:,1:k-1]'*F[:,k]; F[:,k] = F[:,k] - F[:,1:k-1]*dotVF
      end
      # Iterative refinement.
      dotVF = V'*F[:,k]; F[:,k] = F[:,k] - V*dotVF; 
      if k > 1
         dotVF = F[:,1:k-1]'*F[:,k]; F[:,k] = F[:,k] - F[:,1:k-1]*dotVF
      end
      # Orthogonalize F[:,k] against all converged eigenvectors.
      if !isempty(eigvec) 
         F[:,k] = F[:,k] - eigvec*(F[:,k]'*eigvec)'; doteF = (F[:,k]'*eigvec)';
         if norm(doteF) > sqrteps
            F[:,k] = F[:,k] - eigvec*doteF
         end
      end
      F[:,k] = F[:,k]/norm(F[:,k]);
      # Set diagonal element to zero.
      R[k,k] = 0;   
   else
      R[k,k] = norm(F[:,k]);
      F[:,k] = F[:,k] / R[k,k];
   end
   # Modified Gram-Schmidt orthogonalization.
   for j = k+1:blsz
      R[k,j] = dot(F[:,k], F[:,j])
      F[:,j] = F[:,j] - R[k,j] * F[:,k];
   end
   # Iterative refinement
   for j = k+1:blsz
      dotFF = dot(F[:,k], F[:,j])
      R[k,j] = R[k,j] + dotFF;
      F[:,j] = F[:,j] - dotFF * F[:,k];
   end   
end

   return (F,R)

end

#----------------------------------#
# END: BLOCK LANCZOS DECOMPOSITION #
#----------------------------------#

#--------------------------#
# BEGIN: CONVERGENCE TESTS #
#--------------------------#

function convtests(computvec,conv,deflate,dispr,eigval,eigvec,eigresdisp,
                   iter,K,nval,pritz,residuals,ritz,ritzconv,ritzvec,sigma,
                   singblk,sqrteps,tol,Tsz,V)
# This function checks the convergence of Ritz values and Ritz vectors. 

# James Baglama
# DATE: 11/06/01
 
global Rmax;

# Initialization of local variables.
dif  = []; # Place holder for the difference of sets Jr(Jre) and Jrv(Jrev).
Jr   = []; # Place holder for which Ritz values converged.
Jre  = []; # Place holder for which desired Ritz values onverged.
Jrv  = []; # Place holder for which Ritz vectors converged.
Jrev = []; # Place holder for which desired Ritz vectors converged.
ST   = []; # Place holder for which Ritz values stagnated.

# Compute maximum Ritz value to estimate ||A||_2.
if isempty(Rmax)
   Rmax = abs(ritz[Tsz]);
else
   Rmax = max(Rmax,abs(ritz[Tsz]));
end   
Rmax = max(eps()^(2/3),Rmax); 

# Compute tolerance to determine when a Ritz Vector has converged.
RVTol = min(sqrteps,tol); 

# Check for stagnation of Ritz values. eps*100 is used for a tolerance
# to determine when the desired Ritz values are stagnating.
if iter > 1
   if length(pritz) == length(ritz)
      ST = find(abs(pritz-ritz) .< eps()*100);
   end 
end    
       
# Check for convergence of Ritz values and vectors.
   # Smallest eigenvalues.      
   if sigma == "SE"
      
      # Check for convergence of Ritz vectors.
      Jrv  = find(residuals .< RVTol*Rmax);
      Jrev = find(Jrv .<= K-conv);
       
      # Check for convergence of Ritz values.
      Jr = union(find(residuals .< tol*Rmax),ST);
      Jre =  find(Jr  .<= K-conv);
      
      # Output intermediate results.
      if dispr != 0 
         dispeig = [eigval;ritz[1:K-conv]];
         disperr = [eigresdisp';residuals[1:K-conv]'];
         JI = sortperm(dispeig); dispeig = dispeig[JI]; disperr = disperr[JI];
         dispeig = dispeig[1:K]; disperr = disperr[1:K];
         disp(sprintf("      Ritz           Residual       Iteration: %d",iter));
         S = sprintf("%15.5e %15.5e \n",[dispeig';disperr']);
         disp(S); disp(" "); disp(" ");
      end
           
   # Largest eigenvalues.
   elseif sigma == "LE"
             
      # Check for convergence of Ritz vectors.
      Jrv = find(residuals .< RVTol*Rmax);
      Jrev = find(Jrv .>= Tsz-K+1+conv);
            
      # Check for convergence of Ritz values.
      Jr = union(find(residuals .< tol*Rmax),ST);
      Jre =  find(Jr  .>= Tsz-K+1+conv);
      
      
      # Output intermediate results.
      if dispr != 0
         dispeig = [eigval;ritz[Tsz-K+1+conv:Tsz]];
         disperr = [eigresdisp';residuals(Tsz-K+1+conv:Tsz)'];
         JI = sortperm(dispeig); dispeig = dispeig[JI]; disperr = disperr[JI];
         dispeig = dispeig[length(dispeig)-K+1:length(dispeig)]
         disperr = disperr[length(dispeig)-K+1:length(dispeig)]
         #@printf("      Ritz           Residual       Iteration: %d",iter)
         #@printf("%15.5e %15.5e \n",[dispeig';disperr'])
         println()
         println()
      end
       
   # Eigenvalues near nval.
   elseif sigma == "IE"
       
      # Determine a window where the desired Ritz values will occur. 
      JI = sortperm(abs(ritz-nval));
      ritz = ritz[JI]; ritzvec = ritzvec[:,JI];
      residuals = residuals[JI];
     
      # Check for convergence of Ritz vectors.
      Jrv  = find(residuals .< RVTol*Rmax);
      Jrev = find(Jrv .<= K-conv);
       
      # Check for convergence of Ritz values.
      Jr = union(find(residuals .< tol*Rmax),ST);
      Jre =  find(Jr  .<= K-conv);
       
      # Output intermediate results.
      if dispr != 0 
         # Sort output values.
         JI = sortperm(ritz[1:K-conv]);
         sritz = ritz[1:K-conv]; sresiduals =residuals[1:K-conv];
         sritz = sritz[JI]; sresiduals = sresiduals[JI];
         dispeig = [eigval;sritz];
         disperr = [eigresdisp';sresiduals'];
         JI = sortperm(dispeig); dispeig = dispeig[JI]; disperr = disperr[JI];
         dispeig = dispeig[1:K]; disperr = disperr[1:K];
         disp(sprintf("      Ritz           Residual       Iteration: %d",iter));
         S = sprintf("%15.5e %15.5e \n",[dispeig';disperr']);
         disp(S); disp(" "); disp(" ");
      end
end # Switch sigma.

# Remove common values in Jre and Jrev. Common values indicate a desired 
# Ritz pair has converged and will be deflated.
if !isempty(Jr[Jre])
   dif = setdiff(Jr[Jre],Jrv[Jrev]);
end  

# Determine the number of converged desired Ritz vectors.
conv = conv + length(Jrev);

# Determine if the requested number of desired Ritz values have converged.
if conv+length(dif) >= K 
   eigval = [eigval; ritz[union(Jr[Jre],Jrv[Jrev])]];
      
   # If eigenvectors are requested then compute eigenvectors.
   if computvec == "T"
      if isempty(eigvec)
         eigvec = V*ritzvec[:,union(Jr[Jre],Jrv[Jrev])];
      else   
         eigvec = [eigvec V*ritzvec[:,union(Jr[Jre],Jrv[Jrev])]];
      end   
   end   
   
   # Set convergence to true and exit.
   ritzconv = "T"
   return (conv,deflate,eigval,eigvec,eigresdisp,pritz,ritzconv,singblk)
   
end  

# Store previous values of Ritz values to check for stagnation.
I = collect(1:length(ritz))
I = symdiff(I,Jrv)
pritz = ritz[I]

# Compute converged Ritz vectors so they can be deflated.    
if length(Jrv) > 0
   eigval = [eigval;ritz[Jrv]];
   if isempty(eigresdisp)
      eigresdisp = residuals[Jrv]
   else
      new_h = max(size(residuals[Jrv],1), size(eigresdisp,1))
      eigresdisp = [[eigresdisp; zeros(new_h - size(eigresdisp,1), size(eigresdisp,2))] [residuals[Jrv]; zeros(new_h - size(residuals[Jrv],1), 1)]]; 
   end
   if isempty(eigvec)
      eigvec = V*ritzvec[:,Jrv];
   else   
      eigvec = [eigvec V*ritzvec[:,Jrv]];
   end   
   singblk[1] = false;
   
   # Check for convergence of Ritz vectors that do not occur in
   # the desired part of the spectrum. The end points of the 
   # dampening interval(s) will be reset in the main iteraion loop.
   if abs(length(Jrev)-length(Jrv)) != 0
      deflate = true
   end
end

(conv,deflate,eigval,eigvec,eigresdisp,pritz,ritzconv,singblk)
end

#------------------------#
# END: CONVERGENCE TESTS #
#------------------------#

#-----------------------------------------------------------#  
# BEGIN: DETERMINE INTERVALS, LEJA POINTS, AND APPLY SHIFTS #
#-----------------------------------------------------------#

function applyshifts(blsz,endpt,F,fLejapts,iter,K,lcandpts,leftendpt,leftintendpt,
                     lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,rcandpts,
                     rightendpt,rightintendpt,ritz,ritzvec,rprd,rightLejapt,rindex,
                     zertyp,sigma,singblk,sizint,sqrteps,Tsz,T,V,flcount)
# This function determines the endpoints of the dampening intervals, zeros (weighted Leja 
# points or mapped Leja points) and applies the zeros as shifts. 

# James Baglama
# DATE: 8/9/2006 
# (Update the way shifts are applied lines 1254-1269).

# Values used in the GUI demo IRBLDEMO. Not needed for command line computation.
global err, output, waithan, dnbls

# Determine intervals for dampening.
   # Searching for the smallest eigenvalues.
   if sigma == "SE" 
      if isempty(rightendpt)
         rightintendpt  = ritz(Tsz-sizint); rightendpt = ritz(Tsz);
      else
         if !singblk(1) # If a singular block is detected, do not modify end points. 
            if endpt == "MON"
               rightintendpt = min(ritz(Tsz-sizint),rightintendpt);
            else
               rightintendpt  = ritz(Tsz-sizint);
            end   
            rightendpt = max(ritz(Tsz),rightendpt);
         end   
      end
      
      # Choose a value to scale the Leja polynomial to avoid underflow/overflow.
      if isempty(fLejapts)
         norlpol = ritz[1]
      end
      
   # Searching  for the largest eigenvalues.
   elseif sigma == "LE" 
      if isempty(leftendpt) 
         leftendpt = ritz[1];leftintendpt = ritz[1+sizint];
      else
         if !singblk[1] # If a singular block is detected, do not modify end points.
            if endpt == "MON"
               leftintendpt = max(ritz[1+sizint],leftintendpt);
            else
               leftintendpt = ritz[1+sizint];
            end   
            leftendpt = min(ritz[1],leftendpt);
         end  
      end
      
      # Choose a value to scale the Leja polynomial to avoid underflow/overflow.
      if isempty(fLejapts)
         norlpol = ritz[length(ritz)]
      end
      
   # Searching for interior eigenvalues or eignvalues near nval.
   elseif sigma == "IE"
     
      # Shift the spectrum of T by nval.
      T = T - nval*speye(Tsz); ritz = ritz - nval; norlpol = 0;
                            
      # Check the minimum eigenvalue of T.
      (rmin,ri) = findmin(abs(ritz));
      
      # If minimum eigenvalue < sqrt(eps) adjust nval to avoid a singular matrix T.
      if (rmin < sqrteps)
                          
         # If the matrix T is singular adjust nval to be the midpoint between nval 
         # and the next closest Ritz value that it is not equal to. 
         if ((ri > 1) && (ri < Tsz))
            adjust = min(abs(ritz(ri-1)-ritz(ri))/2,abs(ritz(ri+1)-ritz(ri))/2);
         end
         if ri == 1
            adjust = abs(ritz(ri+1)-ritz(ri))/2
         end 
         if ri == Tsz
            adjust = abs(ritz(ri-1)-ritz(ri))/2
         end
         if adjust < sqrteps
            adjust = 2*sqrteps
         end
         T = spdiags(diag(T)-adjust,0,T); ritz = ritz - adjust;
         nval = nval + adjust;
      end   
       
      # Compute the Harmonic Ritz values.  
      B  = sparse(triu(qr(F,0))); B = B[1:blsz,1:blsz];
      EB = sparse(B*rot90(flipud(speye(Tsz,blsz))));
      Hritz = eig(full(T + T\(EB)'*EB)); Hritz = sort(real(Hritz));
       
      # Set extreme endpoints.
      if isempty(leftendpt)
         leftendpt = ritz(1); rightendpt = ritz(Tsz);
      else
         leftendpt = min(leftendpt,ritz(1)); rightendpt = max(rightendpt,ritz(Tsz));
      end   
     
      # Determine interior endpoints.
      kn = 1 + sizint; 
      kp = length(Hritz)-sizint;
           
      while ((kn < length(Hritz)) && (kp > 1))
                  
         # Check to see if the "left" side (negative) interval is acceptable.
         if ((kn < (maximum(find(Hritz.<0)) - round(K/2))) && (Hritz(kn) > leftendpt))
            if endpt == "FLT"
               leftintendpt = Hritz(kn)
            end 
            if endpt == "MON"
               if isempty(leftintendpt)
                  leftintendpt = Hritz(kn); 
               else
                  leftintendpt = max(leftintendpt,Hritz(kn)); 
               end   
            end   
         else
            leftintendpt = []; 
         end
         
         # Check to see if the "right" side (positive) interval is acceptable.
         if ((kp >(minimum(find(Hritz .> 0)) + round(K/2))) && (Hritz(kp) < rightendpt))
            if endpt == "FLT"
               rightintendpt = Hritz(kp)
            end
            if endpt == "MON"
               if isempty(rightintendpt)
                  rightintendpt = Hritz(kp); 
               else
                  rightintendpt = min(rightintendpt,Hritz(kp)); 
               end 
            end   
         else
            rightintendpt = []; 
         end
             
         # Both intervals are unacceptable. Adjust locations and re-compute.
         if ((isempty(leftintendpt)) && (isempty(rightintendpt)))
            kn = kn + 1; kp = kp - 1; 
         else
            break;
         end
         
      end # While loop.
  
      # No interior endpoints, hence dampening intervals, can not be picked. Increase the 
      # number of blocks and restart with BLSZ number of Ritz vectors that are closest
      # to the desired value NVAL. If the matrix T has a singular block, then reset the
      # interior endpoints.
      if ((isempty(leftintendpt)) && (isempty(rightintendpt)) || singblk(1)) 
         lcandpts = []; lindex = []; lprd = []; rcandpts = Float64[]; rindex = []; rprd = Float64[]; fLejapts = Float64[];
         JI = sortperm(abs(ritz));
         if singblk(1) && blsz > 1
            # Compute which vector(s) of V need to be modified.
            k = rem(singblk[2:length(singblk)]-1,blsz)+1; k = union(k,k);
            V[:,k] = randn(size(V,1),length(k));  leftintendpt = []; rightintendpt = [];
         end 
         V = V*ritzvec[:,JI[1:blsz]]; nbls = nbls + 1; dnbls = nbls;
         warning(["Increasing NBLS to ",num2str(nbls)]);
         return (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,
         rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex,V,flcount)
      end   
      
end # Switch sigma.

# Determine shifts. Begin of switch for zertyp.
   # Compute fast Leja points on [-2,2] and map to the dampening interval. This
   # case is not used when sigma = "IE".
   if zertyp == "ML"
      if flcount+nbls > maxdpol
         flcount = 0; maxdpol = length(fLejapts); 
      end
      if flcount == length(fLejapts)
         (fLejapts,rcandpts,rindex,rprd) = fastleja(fLejapts,nbls,rcandpts,rindex,rprd)
      end
      pshifts = fLejapts[flcount+1:nbls+flcount];
      flcount = flcount + length(pshifts);
      if sigma == "SE"
         pshifts = (rightendpt-rightintendpt)/4*(pshifts-2) + rightendpt;
      else
         pshifts = (leftintendpt-leftendpt)/4*(pshifts-2) + leftintendpt;
      end    
                    
   # Compute weighted fast Leja points.     
   elseif zertyp == "WL" 
      if length(fLejapts) + nbls > maxdpol
         # Reset values.
         lcandpts = []; lindex = []; lprd = []; rcandpts = Float64[]; rindex = []; rprd = Float64[]; fLejapts = Float64[];
         leftLejapt = []; rightLejapt = [];
         if sigma == "SE"
            norlpol = ritz(1)
         end
         if sigma == "LE"
            norlpol = ritz(length(ritz))
         end
         if sigma == "IE"
            norlpol = 0
         end
      end 
      (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,pshifts,
      rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex) = 
      fastlejawght(fLejapts,lcandpts,leftendpt,leftintendpt,lprd,
                   leftLejapt,lindex,norlpol,nbls,rcandpts,rightendpt,
                   rightintendpt,rprd,rightLejapt,rindex);
                 
      # Check for early termination and reset values and re-compute fast Leja points.
      if length(pshifts) != nbls
         # Reset values, and try to compute Leja points from the beginning.
         lcandpts = []; lindex = []; lprd = []; rcandpts = Float64[]; rindex = []; rprd = Float64[]; fLejapts = Float64[];
         leftLejapt = []; rightLejapt = [];
         if sigma == "SE"
            norlpol = ritz(1)
         end
         if sigma == "LE"
            norlpol = ritz(length(ritz))
         end
         if sigma == "IE"
            norlpol = 0
         end
         (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,pshifts,
         rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex) = 
         fastlejawght(fLejapts,lcandpts,leftendpt,leftintendpt,lprd,
                      leftLejapt,lindex,norlpol,nbls,rcandpts,rightendpt,
                      rightintendpt,rprd,rightLejapt,rindex);
         # Reset values, if an error has occurred and return using Ritz vectors.
         if length(pshifts) != nbls
            lcandpts = []; lindex = []; lprd = []; rcandpts = Float64[]; rindex = []; rprd = Float64[]; fLejapts = Float64[];
            leftLejapt = []; rightLejapt = [];
            JI = sortperm(abs(ritz));
            V = V*ritzvec[:,JI[1:blsz]];
            return (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,
            rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex,V,flcount)
         end   
      end 
            
end # Switch zertyp.

# Set-up Q to accumulate the Q product.
Q_p = eye(size(T)...);

# Apply nbls-1 shifts.
for J = nbls:-1:2
   (Q,R) = qr(T-pshifts[nbls-J+1]*eye(Tsz));
   T = Q'*(T*Q);
   T=tril(triu(T,-blsz),blsz);
   Q_p = Q_p*Q;
end   
V = V*Q_p;
F = V[:,blsz+1:2*blsz]*T[blsz+1:2*blsz,1:blsz] + F*Q_p[Tsz-blsz+1:Tsz,1:blsz]

# Apply the last shift.  
V = F + V[:,1:blsz]*(T[1:blsz,1:blsz]-pshifts[nbls]*eye(blsz));

# Compute a random vector if a singular block occurs and no eigenvectors converged
# and modify starting vector(s). 
# Details of singular blocks in the IRBL method are given in the paper "Dealing With 
# Linear Dependence during the Iterations of the Restarted Block Lanczos Methods", 
# J. Baglama, Num. Algs., 25, (2000) pp. 23-36. 
# Modify starting vector (case when sigma="SE" and sigma="LE") if a singular block 
# occurs and no eigenvectors converged.  
if singblk[1] != 0 && blsz > 1
   # Compute which vector(s) of V need to be modified.
   k = rem(singblk[2:length(singblk)]-1,blsz)+1; k = union(k,k);
   V[:,k] = randn(size(V,1),length(k)); 
end   

return (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,
         rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex,V,flcount)

end
 
#---------------------------------------------------------#  
# END: DETERMINE INTERVALS, LEJA POINTS, AND APPLY SHIFTS #
#---------------------------------------------------------#

#-------------------------------------------------------------------------------#
# BEGIN: FAST LEJA POINT ROUTINE WITH WEIGHT FUNCTION FOR ONE OR TWO INTERVALS. #
#-------------------------------------------------------------------------------#

function fastlejawght(fLejapts,lcandpts,leftendpt,leftintendpt,lprd,
                      leftLejapt,lindex,norlpol,numshfts,rcandpts,
                      rightendpt,rightintendpt,rprd,rightLejapt,rindex)
# This function computes the fast Leja points with a weight function for one or two intervals. 
#
# References:
# 1.) "Fast Leja Points", J. Baglama, D. Calvetti, and L. Reichel, ETNA (1998) Volume 7.
# 2.) "Iterative Methods for the Computation of a Few Eigenvalues of a Large Symmetric Matrix",
#      J. Baglama, D. Calvetti, and L. Reichel, BIT, 36 (1996), pp. 400-421.

# James Baglama 
# DATE: 11/06/01

# Values used in the GUI demo IRBLDEMO. Not needed for command line computation.
global err, output, waithan

# Number of previously generated fast Leja points.
nfLejapts = length(fLejapts); 
    
# Initialize arrays.
pshifts = []; # Output array that holds the shifts to be applied. 
IL = [];      # Place holder for finding points in the left intervals.
IR = [];      # Place holder for finding points in the right intervals.

# Checking for errors in input argument(s). 
if (numshfts < 1) 
   err = "Number of shifts is (<= 0) in fast Leja with weight function.";
   error(err)
end
if ((isempty(leftintendpt)) && (isempty(rightintendpt)))
   err = "Empty endpoints in fast Leja with weight function.";
   error(err)
end

# Initialize the arrays fLejapts, lcandpts, and rcandpts.
if ((isempty(rcandpts)) && (!isempty(rightintendpt)))
   nfLejapts = nfLejapts+1;
   fLejapts[nfLejapts] = rightendpt; rightLejapt = nfLejapts; 
   rcandpts[1] = (rightintendpt+fLejapts(nfLejapts))/2;
   rprd[1] = prod((rcandpts[1]-fLejapts)./(norlpol-fLejapts));
   rindex[1,1] = 0; rindex[1,2] = nfLejapts;
   numshfts = numshfts-1; pshifts = fLejapts[nfLejapts];
end   
if ((isempty(lcandpts)) && (!isempty(leftintendpt)))
   nfLejapts = nfLejapts+1;
   fLejapts[nfLejapts] = leftendpt; leftLejapt = nfLejapts;   
   lcandpts[1] = (leftintendpt+fLejapts(nfLejapts))/2;     
   lprd[1] = prod((lcandpts[1]-fLejapts)./(norlpol-fLejapts));
   lindex[1,1] = nfLejapts; lindex[1,2] = 0;
   numshfts = numshfts-1; pshifts = [pshifts fLejapts[nfLejapts]];
end

# Set products of the end points to zero.
rightprd = 0; leftprd = 0;

# Find closest Leja points in each interval to the interior end points.
sizlcpts = length(lcandpts); sizrcpts = length(rcandpts);    
if ((sizlcpts > 1) && (!isempty(leftintendpt)))
   IL = find(fLejapts .< leftintendpt);(YL,KL) = findmax(fLejapts[IL]);      
   if !isempty(IL(KL))
      lcandpts[1] = (leftintendpt+fLejapts(IL(KL)))/2;     
      lprd[1] = prod((lcandpts[1]-fLejapts)./(norlpol-fLejapts));
      lindex[1,1] = IL(KL); lindex[1,2] = 0;
   else
      lprd[1] == 0;
   end
   leftprd = prod((leftendpt-fLejapts)./(norlpol-fLejapts)); 
end
if ((sizrcpts > 1) && (!isempty(rightintendpt)))
   IR = find(fLejapts .> rightintendpt);(YR,KR) = findmin(fLejapts[IR]);      
   if !isempty(IR(KR))
      rcandpts[1] = (rightintendpt+fLejapts(IR(KR)))/2;     
      rprd[1] = prod((rcandpts[1]-fLejapts)./(norlpol-fLejapts));
      rindex[1,1] = 0; rindex[1,2] = IR(KR);
   else
      rprd[1] = 0;
   end
   rightprd = prod((rightendpt-fLejapts)./(norlpol-fLejapts));
end

# Initialize values.
JL = [];        # Index of acceptable candidate points in the left interval.
JR = [];        # Index of acceptable candidate points in the right interval.
maxprdJL = -1;  # Maximum value of polynomial on left interval.
maxprdJR = -1;  # Maximum value of polynomial on right interval.
maxJL = -1;     # Index of maximum value of polynomial on left interval.
maxJR = -1;     # Index of maximum value of polynomial on right interval.
  
# Search for Leja points. 
for k = 1+nfLejapts:numshfts+nfLejapts
   
   if !isempty(leftintendpt)
      JL = find(lcandpts .> leftendpt && lcandpts < leftintendpt);
      (maxprdJL,maxJL) = findmax(abs((lcandpts(JL)-leftintendpt).*lprd(JL)));
      maxJL = JL(maxJL); leftprd = abs((leftendpt-leftintendpt).*leftprd);
      if isempty(maxprdJL)
         maxprdJL = leftprd; maxJL = 0
      elseif(leftprd > maxprdJL)
         maxprdJL = leftprd; maxJL = 0
      end
   end
      
   if !isempty(rightintendpt)
      JR = find(rcandpts .> rightintendpt && rcandpts < rightendpt);
      (maxprdJR,maxJR) = findmax(abs((rcandpts[JR]-rightintendpt).*rprd(JR)));
      maxJR = JR(maxJR); rightprd = abs((rightendpt-rightintendpt).*rightprd);
      if isempty(maxprdJR)
         maxprdJR = rightprd; maxJR = 0
      elseif rightprd > maxprdJR
         maxprdJR = rightprd; maxJR = 0
      end
   end
   
   # Check for overflow or underflow.
   if ((isinf([maxprdJR,maxprdJL])) || (isnan([maxprdJR,maxprdJL])))
      lcandpts = []; rcandpts = Float64[]; fLejapts = Float64[]; lindex = []; rindex = []; lprd = []; rprd = Float64[];
      rightLejapt = []; leftLejapt = []; pshifts = []; return; 
   end    
    
   # Case 1. Right end point is the next fast Leja point.
   if ((maxprdJR >= maxprdJL) && (maxJR == 0))
      fLejapts[k] = rightendpt; sizrcpts = sizrcpts+1;
      if isempty(IR)
         rcandpts[sizrcpts] = (rightendpt + rightintendpt)/2;
         rindex[sizrcpts,1] = 0; rindex[sizrcpts,2] = k;
      else   
         rcandpts[sizrcpts] = (rightendpt + fLejapts[rightLejapt])/2;
         rindex[sizrcpts,1] = rightLejapt; rindex[sizrcpts,2] = k;
      end   
      rightLejapt = k; rightprd = 0; 
      rprd[sizrcpts] = prod((rcandpts[sizrcpts]-fLejapts[1:k-1])./(norlpol-fLejapts[1:k-1]));
   end   
   
   # Case 2. Left end point is the next fast Leja point.
   if ((maxprdJL > maxprdJR) && (maxJL == 0))
      fLejapts[k] = leftendpt; sizlcpts = sizlcpts+1;
      if isempty(IL)
         lcandpts(sizlcpts) = (leftendpt + leftintendpt)/2;
         lindex[sizlcpts,1] = k; lindex[sizlcpts,2] = 0;
      else   
         lcandpts(sizlcpts) = (leftendpt + fLejapts[leftLejapt])/2;
         lindex[sizlcpts,1] = k; lindex[sizlcpts,2] = leftLejapt;
      end
      leftLejapt = k; leftprd = 0;
      lprd(sizlcpts) = prod((lcandpts(sizlcpts)-fLejapts[1:k-1])./(norlpol-fLejapts[1:k-1]));
   end   
   
   # Case 3. Leja point picked in Left interval.
   if (maxprdJL > maxprdJR) && (maxJL != 0)
      fLejapts[k] = lcandpts(maxJL); sizlcpts = sizlcpts+1; 
      lindex[sizlcpts,1] = k; lindex[sizlcpts,2] = lindex[maxJL,2]; lindex[maxJL,2] = k;
      if (lindex[sizlcpts,2] == 0)
         lcandpts(sizlcpts) = (fLejapts[lindex[sizlcpts,1]]+leftintendpt)/2;
      else 
         lcandpts(sizlcpts) = (fLejapts[lindex[sizlcpts,1]]+fLejapts[lindex[sizlcpts,2]])/2;
      end
      lcandpts(maxJL) = (fLejapts[lindex[maxJL,1]]+fLejapts[lindex[maxJL,2]])/2;
      lprd(maxJL) = prod((lcandpts(maxJL)-fLejapts[1:k-1])./(norlpol-fLejapts[1:k-1]));
      lprd(sizlcpts) = prod((lcandpts(sizlcpts)-fLejapts[1:k-1])./(norlpol-fLejapts[1:k-1]));
   end       
      
   # Case 4. Leja point picked in right interval. 
   if ((maxprdJR >= maxprdJL) && (maxJR != 0))
      fLejapts[k] = rcandpts[maxJR]; sizrcpts = sizrcpts+1; 
      rindex[sizrcpts,1] = k; rindex[sizrcpts,2] = rindex[maxJR,2]; rindex[maxJR,2] = k;
      if (rindex[maxJR,1] == 0)
         rcandpts[maxJR] = (rightintendpt+fLejapts(rindex[maxJR,2]))/2;
      else 
         rcandpts[maxJR] = (fLejapts[rindex[maxJR,1]]+fLejapts[rindex[maxJR,2]])/2;
      end
      rcandpts[sizrcpts] = (fLejapts[rindex[sizrcpts,1]]+fLejapts[rindex[sizrcpts,2]])/2;
      rprd[maxJR] = prod((rcandpts[maxJR]-fLejapts[1:k-1])./(norlpol-fLejapts[1:k-1]));
      rprd[sizrcpts] = prod((rcandpts[sizrcpts]-fLejapts[1:k-1])./(norlpol-fLejapts[1:k-1]));
   end 
   rprd = rprd.*(rcandpts-fLejapts(k))/(norlpol-fLejapts(k));
   lprd = lprd.*(lcandpts-fLejapts(k))/(norlpol-fLejapts(k));
   if rightprd != 0
      rightprd = rightprd.*(rightendpt - fLejapts[k])/(norlpol-fLejapts(k))
   end
   if leftprd != 0
      leftprd =  leftprd.*(leftendpt - fLejapts[k])/(norlpol-fLejapts(k))
   end

end # For loop.

pshifts = [pshifts,fLejapts[nfLejapts+1:nfLejapts+numshfts]];

return (fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,pshifts,
         rcandpts,rightendpt,rightintendpt,rprd,rightLejapt,rindex)
end

#-----------------------------------------------------#
# END: FAST LEJA POINT ROUTINE WITH WEIGHT FUNCTION #
#-----------------------------------------------------#

#-------------------------------------------#
# BEGIN: FAST LEJA POINT ROUTINE FOR [-2,2] #
#-------------------------------------------#

function fastleja(fLejapts,numshfts,rcandpts,rindex,rprd)
# This function is made to be used in irbleigs.m. This function computes the
# fast Leja points on the interval [-2,2]. 
#
# References:
# "Fast Leja Points", J. Baglama, D. Calvetti, and L. Reichel, ETNA (1998) Volume 7.

# James Baglama
# 10/18/01

# Values used in the GUI demo IRBLDEMO. Not needed for command line computation. 
global err, output, waithan

# Number of previously generated fast Leja points.
nfLejapts = length(fLejapts)
    
# Checking for errors in input argument(s). 
if (numshfts < 1) 
   err = "Number of shifts is 0 or (< 0) in fast Leja routine.";
   error(err)
end
   
# Initialize the array fLejapts.
if (nfLejapts < 2)
   push!(fLejapts, 2)
   if ((nfLejapts == 0) && (numshfts == 1))
      return
   end 
   push!(fLejapts, -2)
   push!(rcandpts, (fLejapts[1]+fLejapts[2])/2)
   push!(rprd, prod(rcandpts[1]-fLejapts))
   rindex = [ 2 1 ]
   if ((nfLejapts == 0) && (numshfts == 2))
      return
   end
   if ((nfLejapts == 1) && (numshfts == 1))
      return
   end
   nfLejapts = 2; numshfts = numshfts-2;
end    
    
# Main loop to calculate the fast Leja points.
for k = 1+nfLejapts:numshfts+nfLejapts
   (maxval,maxk)  = findmax(abs(rprd));   
   push!(fLejapts, rcandpts[maxk]);
   rindex = [rindex; k rindex[maxk,2]]
   rindex[maxk,2] = k;
   rcandpts[maxk] = (fLejapts[rindex[maxk,1]]+fLejapts[rindex[maxk,2]])/2;
   push!(rcandpts, (fLejapts[rindex[k-1,1]]+fLejapts[rindex[k-1,2]])/2)
   rprd[maxk]     = prod(rcandpts[maxk]-fLejapts[1:k-1]);
   push!(rprd, prod(rcandpts[k-1]-fLejapts[1:k-1]))
   rprd           = rprd.*(rcandpts-fLejapts[k]);
end

return (fLejapts,rcandpts,rindex,rprd)
end

#--------------------------------------------#
# END: FAST LEJA POINT ROUTINE ONE INTERVAL. #
#--------------------------------------------#

srand(1234)
A=matrix_market_read(ARGS[1], true, true)

eigs(A;nev=5)

eigs_time = -time()
(d,nconv,niter,nmult,resid) = eigs(A;nev=5)
eigs_time += time()
println(d)

irbleigs(A)

global mainloop_time = 0.
global deflate_time = 0.
global blanz_time = 0.
global eig_time = 0.
global sort_time = 0.
global residual_time = 0.
global convtests_time = 0.
global applyshifts_time = 0.

global blanz_ortho_time = 0.
global blanz_spmv_time = 0.
global blanz_ortho_inner_time = 0.
global blanz_qr_time = 0.
global blanz_etc1_time = 0.
global blanz_qrsingblk_time = 0.
global blanz_etc2_time = 0.

global dgemm1_time = 0.
global dgemm2_time = 0.
global dgemm3_time = 0.
global dgemm4_time = 0.
global dgemm5_time = 0.

irbleigs_time = -time()
d2, iter = irbleigs(A)
irbleigs_time += time()
println(d2)

println("iter = $iter")
println("main_loop = $mainloop_time deflate = $deflate_time blanz = $blanz_time eig = $eig_time sort = $sort_time residual = $residual_time convtests = $convtests_time applyshifts = $applyshifts_time")
println("blanz_ortho_time = $blanz_ortho_time blanz_spmv_time = $blanz_spmv_time blanz_ortho_inner_time = $blanz_ortho_inner_time blanz_qr_time = $blanz_qr_time blanz_etc1_time = $blanz_etc1_time blanz_qrsingblk_time = $blanz_qrsingblk_time blanz_etc2_time = $blanz_etc2_time")
println("dgemm1_time = $dgemm1_time dgemm2_time = $dgemm2_time dgemm3_time = $dgemm3_time dgemm4_time = $dgemm4_time dgemm5_time = $dgemm5_time")

println("eigs = $eigs_time irbleigs = $irbleigs_time")

# vim: set tabstop=8 softtabstop=3 sw=3 expandtab:
