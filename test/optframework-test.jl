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

call_sig_arg_tuple = (Array{Float64,1}, Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, Array{Float64,1}, Float64, Int64)
ast = code_lowered(cg, call_sig_arg_tuple)
println("******************* lowered AST for cg **************")
println(ast)

# code copied and modified from OptFramework.jl
method = Base._methods(cg, call_sig_arg_tuple, -1)
assert(length(method) == 1)
method = method[1]

println("Initial code to optimize = ", ast)

method[3].isstaged = true
method[3].func.code.ast = ccall(:jl_compress_ast, Any, (Any,Any), method[3].func.code, ast)
linfo = Base.func_for_method(method[3], call_sig_arg_tuple, method[2])
# Must be going from lowered AST to type AST.
(cur_ast, ty) = Base.typeinf(linfo, method[1], method[2])

