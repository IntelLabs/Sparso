include("utils.jl")

function pcg(x, A, b, M, tol, maxiter)
    r = b - A * x
    normr0 = norm(r)
    rel_err = 1.0
    z::Vector{Float64} = M \ r   # Somehow, types for "M\r" is not inferred. 
                                 # The explicit typing makes it compile and shows
                                 # the call replacement, but has a running error 
                                 # in type conversion. So this example is for 
                                 # manual checking of call replacement only.
    p = copy(z) #NOTE: do not write "p=z"! That would make p and z aliased (the same variable)
    rz = dot(r, z)
    k = 1
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
        z = M \ r  
        rz = dot(r, z)
        beta = rz/old_rz
        p = z + beta * p
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
M = A # Perfect
