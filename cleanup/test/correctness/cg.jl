include("utils.jl")

function cg(x, A, b, tol, maxiter)
    r = b - A * x
    rel_err = 1
    p = copy(r) #NOTE: do not write "p=r"! That would make p and r aliased (the same variable)
    rz = dot(r, r)
    normr0 = sqrt(rz)
    k = 1
    while k <= maxiter
        old_rz = rz
        Ap = A*p #Ap = SparseAccelerator.SpMV(A, p) # manual # This takes most time. Compiler can reorder A to make faster
        alpha = old_rz / dot(p, Ap) #alpha = old_rz / SparseAccelerator.Dot(p, Ap) # manual
        x += alpha * p #SparseAccelerator.WAXPBY!(x, alpha, p, 1, x) # manual
        r -= alpha * Ap #SparseAccelerator.WAXPBY!(r, -alpha, Ap, 1, r) # manual
        rz = dot(r, r) #rz = SparseAccelerator.Dot(r, r) # manual
        rel_err = sqrt(rz)/normr0
        if rel_err < tol 
            break
        end
        beta = rz/old_rz
        p = r + beta * p #SparseAccelerator.WAXPBY!(p, 1, r, beta, p) # manual
        k += 1
    end
    return x, k, rel_err
end

m = 10
A = generate_symmetric_sparse_matrix(m)
x = repmat([1/m], m)
b   = ones(Float64, m)
tol = 1e-10
maxiter = 1000