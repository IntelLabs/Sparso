function pcg_jacobi(x, A, b, tol, maxiter)
    inv_d = 1./diag(A)
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1
    z = inv_d .* r
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
    while k <= maxiter
        old_rz = rz
        Ap = A*p #Ap = SparseAccelerator.SpMV(A, p) # manual # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap) #alpha = old_rz / SparseAccelerator.dot(p, Ap) # manual
        x += alpha * p #SparseAccelerator.WAXPBY!(x, alpha, p, 1, x) # manual
        r -= alpha * Ap #SparseAccelerator.WAXPBY!(r, -alpha, Ap, 1, r) # manual
        rel_err = sqrt(dot(r, r))/normr0 #rel_err = sqrt(SparseAccelerator.dot(r, r))/normr0 # manual
        if rel_err < tol 
            break
        end
        z = inv_d .* r #z = SparseAccelerator.element_wise_multiply(inv_d, r)
        rz = dot(r, z) #rz = SparseAccelerator.dot(r, z) # manual
        beta = rz/old_rz
        p = z + beta * p #SparseAccelerator.WAXPBY!(p, 1, z, beta, p) # manual
        k += 1
    end
    return x, k, rel_err
end