include("../../src/SparseAccelerator.jl")
using SparseAccelerator

module X
    module Y
        module Z
            module U
                module V
                    module W
                        function f(x)
                            x + 1
                        end
                    end
                end
            end
        end
    end
end

function foo()
    X.Y.Z.U.V.W.f(1)
end

function bar(A, B)
    A * B
end

ast = code_typed(foo, (), optimize = false)
call_args = ast[1].args[3].args[2].args[1].args
println("Call args: ", call_args)
module_name, function_name = SparseAccelerator.resolve_call_names(call_args)
println("Module name: ", module_name)
println("Function name: ", function_name)

ast = code_typed(bar, (SparseMatrixCSC, Vector), optimize = false)
call_args = ast[1].args[3].args[2].args[1].args
println("\nCall args: ", call_args)
module_name, function_name = SparseAccelerator.resolve_call_names(call_args)
println("Module name: ", module_name)
println("Function name: ", function_name)
