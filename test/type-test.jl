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


function pcg(x, A, b, M, tol, maxiter)
    r = b - A * x
    z = M \ r
    p = z
    rz = dot(r, z)
    k = 1
    while k <= maxiter
        old_rz = rz
        α = old_rz / dot(p, A * p)
        x += α * p
        r -= α * A * p 
        if norm(r) < tol 
            break
        end
        z = M \ r  
        rz = dot(r, z)
        β = rz/old_rz
        p = z + β * p
        k += 1
    end
    return x #, k
end


ast = code_typed(pcg, (Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Float64, Int64), optimize=false)
println("******************* typed AST for pcg *************")
println(ast)

ast = code_typed(cg, (Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Array{Float64,1}, Float64, Int64), optimize=false)
println("******************* typed AST for cg **************")
println(ast)

