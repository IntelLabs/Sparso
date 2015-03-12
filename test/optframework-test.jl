function cg(x, A)
    A * x
end

call_sig_arg_tuple = (Int64, Int64)
ct = code_typed(cg, call_sig_arg_tuple)
ct = ct[1]
println("******************* typed AST for cg **************")
println(ct)
ast = code_lowered(cg, call_sig_arg_tuple)
ast = ast[1]
println("type of ast = ", typeof(ast))
println("******************* lowered AST for cg **************")
println(ast)

# code copied and modified from OptFramework.jl
method = Base._methods(cg, call_sig_arg_tuple, -1)
assert(length(method) == 1)
method = method[1]

println("Initial code to optimize = \n", ast)

println("second = ", ast.args[1]..., " type = ", typeof(ast.args[1]))
#println("second = ", ast.args[1]..., " ", eval(Expr(:tuple, ast.args[1]...)))
new_func = Expr(:function, Expr(:call, :new_func_name, ast.args[1]...), Expr(:block, ast.args[3].args...))
#new_func = Expr(:function, Expr(:tuple, ast.args[1]...), Expr(:block, ast.args[3].args...))
eval_new_func = eval(new_func)
println("eval_new_func = ", eval_new_func, " type = ", typeof(eval_new_func))

cur_ast = code_typed(eval_new_func, call_sig_arg_tuple)
println("cur_ast = ", cur_ast, " type = ", typeof(cur_ast))
cur_ast = cur_ast[1]

if false
#method[3].func.code.ast = ast
#println("head = ", ast.head)
#println("args3 = ", ast.args[3], " type = ", typeof(ast.args[3]))
stmt = ast.args[3].args[2]
println("stmt = ", stmt)
stmt = stmt.args[1]
println("stmt = ", stmt, " head = ", stmt.head)
arg1 = stmt.args[2]
arg2 = stmt.args[3]
println("arg1 = ", arg1, " type = ", typeof(arg1))
println("arg2 = ", arg2, " type = ", typeof(arg2))

method[3].isstaged = true
#method[3].func.code.ast = ccall(:jl_compress_ast, Any, (Any,Any), method[3].func.code, ast)
linfo = Base.func_for_method(method[3], call_sig_arg_tuple, method[2])
# Must be going from lowered AST to type AST.
(cur_ast, ty) = Base.typeinf(linfo, method[1], method[2])

if !isa(cur_ast,Expr)
    cur_ast = ccall(:jl_uncompress_ast, Any, (Any,Any), linfo, cur_ast)
end

end

println("Typed code after = ", cur_ast, " type = ", typeof(cur_ast))
