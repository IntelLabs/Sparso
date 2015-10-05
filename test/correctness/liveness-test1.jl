include("../../src/SparseAccelerator.jl")
using SparseAccelerator
using CompilerTools
using CompilerTools.CFGs
using CompilerTools.OptFramework
using CompilerTools.LivenessAnalysis

function dump_liveness(func_ast :: Expr, func_arg_types :: Tuple, func_args)
    assert(typeof(func_ast)== Expr)
    assert(func_ast.head == :lambda)

    LivenessAnalysis.set_use_inplace_naming_convention()
    liveness = LivenessAnalysis.from_expr(func_ast)#, no_mod = SparseAccelerator.create_unmodified_args_dict())
    println("Liveness:\n", liveness)

    func_ast
end

sparse_pass = OptFramework.OptPass(dump_liveness, true)
OptFramework.setOptPasses([sparse_pass])
#SparseAccelerator.set_debug_level(2)

include("./cg.jl")

A       = matrix_market_read(ARGS[1], true, true)
m       = size(A, 1)
x       = zeros(Float64, m)
b       = ones(Float64, m)
tol     = 1e-10
maxiter = 1000
@acc x, k, rel_err = cg(x, A, b, tol, maxiter)
