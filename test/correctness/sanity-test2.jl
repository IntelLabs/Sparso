# Testing OptFramework

using CompilerTools.OptFramework

function do_nothing(ast, call_sig_arg_tuple, call_sig_args)
    println("\tI do nothing!")
    ast
end

function foo(A, x)
    A*x
end

pass = OptFramework.OptPass(do_nothing, true)
OptFramework.setOptPasses([pass])

m = 10
A = sprand(m, m, 0.1)
x = repmat([1/m], m)

println("Testing OptFramework:")
@acc foo(A, x)
