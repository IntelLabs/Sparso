function pcg_symgs(x, A, b, tol, maxiter)
    L = tril(A)
    U  = spdiagm(1./diag(A))*triu(A)
    M = L*U
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1

    z = copy(r)
    Base.SparseMatrix.fwdTriSolve!(L, z)
    Base.SparseMatrix.bwdTriSolve!(U, z)

    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    time1 = time()
    while k <= maxiter
        old_rz = rz
        Ap = A*p # Ap = SparseAccelerator.SpMV(A, p) # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        if rel_err < tol 
            break
        end

        z = copy(r)
        Base.SparseMatrix.fwdTriSolve!(L, z)
          # Could have written as z=L\z if \ is specialized for triangular
        Base.SparseMatrix.bwdTriSolve!(U, z)
          # Could have wrriten as z=U\z if \ is specialized for triangular

        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
        k += 1
    end
    return x, k, rel_err
end