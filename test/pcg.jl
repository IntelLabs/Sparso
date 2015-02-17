using MatrixMarket
using Base.Test

function cg(x, A, b, tol, maxiter)
    r = b - A * x
    rel_err = 1
    p = r
    rz = dot(r, r)
    normr0 = sqrt(rz)
    k = 1
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
    return x, k, rel_err
end

function pcg_jacobi(x, A, b, tol, maxiter)
    inv_d = 1./diag(A)

    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = inv_d .* r
    p = z
    rz = dot(r, z)
    k = 1
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

    p = z
    rz = dot(r, z)
    k = 1
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
    return x, k, rel_err
end

function pcg(x, A, b, M, tol, maxiter)
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = M \ r
    p = z
    rz = dot(r, z)
    k = 1
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
    return x, k, rel_err
end

matrices = [
#"hpcg_4",
"hpcg_32",
]

tol = 1e-6
maxiter = 1000

for matrix in matrices
  println(matrix)
  A = MatrixMarket.mmread("data/$matrix.mtx")

  N   = size(A, 1)
  b   = ones(Float64, N)

  x   = zeros(Float64, N)
  @time x, k, err = cg(x, A, b, tol, maxiter)
  # The line below invokes generic PCG.
  # Ideally, want to get same performance with specialized ver. in above line
  #@time x, k, err = pcg(x, A, b, speye(N), tol, maxiter)
  println("Identity preconditioner: $k iterations $err error")

  x   = zeros(Float64, N)
  @time x, k, err = pcg_jacobi(x, A, b, tol, maxiter)
  #@time x, k, err = pcg(x, A, b, spdiagm(diag(A)), tol, maxiter)
  println("Jacobi preconditioner: $k iterations $err error")

  x   = zeros(Float64, N)
  @time x, k, err = pcg_symgs(x, A, b, tol, maxiter)
  #@time x, k, err = pcg(x, A, b, tril(A)*spdiagm(1./diag(A))*triu(A), tol, maxiter)
  println("SymGS preconditioner: $k iterations $err error")

  M = A # Perfect
  x   = zeros(Float64, N)
  @time x, k, err = pcg(x, A, b, M, tol, maxiter)
  println("Perfect preconditioner: $k iterations $err error")
  println()
end
