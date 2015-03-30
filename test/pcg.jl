#TODO: remove this include once the package is ready
include("../src/OptFramework.jl")
include("../src/SparseAccelerator.jl")
include("../src/sparse-analyze.jl")

using OptFramework
using MatrixMarket
using Base.Test

sparse_pass = OptFramework.optPass(SparseAccelerator.SparseOptimize, true)
OptFramework.setOptPasses([sparse_pass])

function cg(x, A, b, tol, maxiter)
    tic()
    r = b - A * x
    rel_err = 1
    p = copy(r) #NOTE: do not write "p=r"! That would make p and r aliased (the same variable)
    rz = dot(r, r)
    normr0 = sqrt(rz)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rz = dot(r, r)
        rel_err = sqrt(rz)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end
        beta = rz/old_rz
        p = r + beta * p
        k += 1
    end
    time2 = time()
    println("Time of original loop= ", time2 - time1, " seconds")
    toc()
    return x, k, rel_err
end

function pcg_jacobi(x, A, b, tol, maxiter)
    tic()
    inv_d = 1./diag(A)

    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = inv_d .* r
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end
        z = inv_d .* r  
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    time2 = time()
    println("Time of original loop= ", time2 - time1, " seconds")
    toc()
    return x, k, rel_err
end

#Check that matrix is square
function chksquare(A::AbstractMatrix)
    m,n = size(A)
    m == n || throw(DimensionMismatch("matrix is not square"))
    m
end

function fwdTriSolve!(A::SparseMatrixCSC, B::AbstractVecOrMat)
# forward substitution for CSC matrices
    n = length(B)
    if isa(B, Vector)
        nrowB = n
        ncolB = 1
    else
        nrowB, ncolB = size(B)
    end
    ncol = chksquare(A)
    if nrowB != ncol
        throw(DimensionMismatch("A is $(ncol)X$(ncol) and B has length $(n)"))
    end

    aa = A.nzval
    ja = A.rowval
    ia = A.colptr

    joff = 0
    for k = 1:ncolB
        for j = 1:(nrowB-1)
            jb = joff + j
            i1 = ia[j]
            i2 = ia[j+1]-1
            B[jb] /= aa[i1]
            bj = B[jb]
            for i = i1+1:i2
                B[joff+ja[i]] -= bj*aa[i]
            end
        end
        joff += nrowB
        B[joff] /= aa[end]
    end
    return B
end

function bwdTriSolve!(A::SparseMatrixCSC, B::AbstractVecOrMat)
# backward substitution for CSC matrices
    n = length(B)
    if isa(B, Vector)
        nrowB = n
        ncolB = 1
    else
        nrowB, ncolB = size(B)
    end
    ncol = chksquare(A)
    if nrowB != ncol throw(DimensionMismatch("A is $(ncol)X$(ncol) and B has length $(n)")) end

    aa = A.nzval
    ja = A.rowval
    ia = A.colptr

    joff = 0
    for k = 1:ncolB
        for j = nrowB:-1:2
            jb = joff + j
            i1 = ia[j]
            i2 = ia[j+1]-1
            B[jb] /= aa[i2]
            bj = B[jb]
            for i = i2-1:-1:i1
                B[joff+ja[i]] -= bj*aa[i]
            end
        end
        B[joff+1] /= aa[1]
        joff += nrowB
    end
   return B
end

function pcg_symgs(x, A, b, tol, maxiter)
    tic()
    L = tril(A)
    U = spdiagm(1./diag(A))*triu(A)
    M = L*U
    #inv_d = 1./diag(A)

    r = b - A * x
    normr0 = norm(r)
    rel_err = 1

    z = copy(r)
    fwdTriSolve!(L, z)
    bwdTriSolve!(U, z)
    #z = M \ r

    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end

        z = copy(r)
        fwdTriSolve!(L, z)
          # Could have written as z=L\z if \ is specialized for triangular
        bwdTriSolve!(U, z)
          # Could have wrriten as z=U\z if \ is specialized for triangular

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    time2 = time()
    println("Time of original loop= ", time2 - time1, " seconds")
    toc()
    return x, k, rel_err
end

function pcg(x, A, b, M, tol, maxiter)
    tic()
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = M \ r
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        #println(rel_err)
        if rel_err < tol 
            break
        end
        z = M \ r  
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    time2 = time()
    println("Time of original loop= ", time2 - time1, " seconds")
    toc()
    return x, k, rel_err
end

#matrices = [
#"hpcg_4",
#"hpcg_32",
#]

tol = 1e-6
maxiter = 1000

#for matrix in matrices
#  println(matrix)
#  A = MatrixMarket.mmread("data/$matrix.mtx")

  A = MatrixMarket.mmread(ARGS[1])
#  A = spones(A)

  N   = size(A, 1)
  b   = ones(Float64, N)

  tests = 5
  
  println("\n\n**** original cg perf")
  for i = 1:tests
      x   = zeros(Float64, N)
      x, k, err = cg(x, A, b, tol, maxiter)
      # The line below invokes generic PCG.
      # Ideally, want to get same performance with specialized ver. in above line
      #@time x, k, err = pcg(x, A, b, speye(N), tol, maxiter)
      println("Identity preconditioner: $k iterations $err error")
  end

  println("\n\n**** accelerated cg perf")
  for i = 1:tests
      x   = zeros(Float64, N)
      @acc result = cg(x, A, b, tol, maxiter)
      x, k, err = result
      # The line below invokes generic PCG.
      # Ideally, want to get same performance with specialized ver. in above line
      #@time x, k, err = pcg(x, A, b, speye(N), tol, maxiter)
      println("Identity preconditioner: $k iterations $err error")
  end

  println("\n\n**** original pcg_jacobi perf")
  for i = 1:tests
      x   = zeros(Float64, N)
      x, k, err = pcg_jacobi(x, A, b, tol, maxiter)
      #@time x, k, err = pcg(x, A, b, spdiagm(diag(A)), tol, maxiter)
      println("Jacobi preconditioner: $k iterations $err error")
  end

  println("\n\n**** accelerated pcg_jacobi perf")
  for i = 1:tests
      x   = zeros(Float64, N)
      @acc result = pcg_jacobi(x, A, b, tol, maxiter)
      x, k, err = result
      #@time x, k, err = pcg(x, A, b, spdiagm(diag(A)), tol, maxiter)
      println("Jacobi preconditioner: $k iterations $err error")
  end
  
  println("\n\n**** original pcg_symgs perf")
  for i = 1:tests
      x   = zeros(Float64, N)
      x, k, err = pcg_symgs(x, A, b, tol, maxiter)
      #@time x, k, err = pcg(x, A, b, tril(A)*spdiagm(1./diag(A))*triu(A), tol, maxiter)
      println("SymGS preconditioner: $k iterations $err error")
  end

  println("\n\n**** accelerated pcg_symgs perf")
  for i = 1:tests
#      x   = zeros(Float64, N)
#      @acc result = pcg_symgs(x, A, b, tol, maxiter)
#      x, k, err = result
      #@time x, k, err = pcg(x, A, b, tril(A)*spdiagm(1./diag(A))*triu(A), tol, maxiter)
#      println("SymGS preconditioner: $k iterations $err error")
  end
  
  println("\n\n**** original pcg perf")
  for i = 1:tests  
      M = A # Perfect
      x   = zeros(Float64, N)
      x, k, err = pcg(x, A, b, M, tol, maxiter)
      println("Perfect preconditioner: $k iterations $err error")
  end

  println("\n\n**** accelerated pcg perf")
  for i = 1:tests  
      M = A # Perfect
      x   = zeros(Float64, N)
      @acc result = pcg(x, A, b, M, tol, maxiter)
      x, k, err = result
      println("Perfect preconditioner: $k iterations $err error")
  end
  
  println()
#end
