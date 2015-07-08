using CompilerTools
using CompilerTools.OptFramework

include("../Src/SparseAccelerator.jl")
using SparseAccelerator

function checkDistributivity(ast, call_sig_arg_tuple, call_sig_args)
  assert(typeof(ast)== Expr)
  assert(ast.head == :lambda)

  symbolInfo = SparseAccelerator.initSymbol2TypeDict(ast)
  distributive = SparseAccelerator.checkDistributivity(ast, symbolInfo, true)
  ast
end

sparse_pass = OptFramework.optPass(checkDistributivity, true)
OptFramework.setOptPasses([sparse_pass])
SparseAccelerator.set_debug_level(2)

function test(A, x)
    v = SparseAccelerator.SpMV(A, x)
    v = A * x
end

m = 10
A = sprand(m, m, 0.1)
x = repmat([1/m], m)
@acc test(A, x)
