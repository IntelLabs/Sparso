module StructureAnalysisStructOnly
   
    using SparseAccelerator: Sym, Symexpr, TypedExprNode
    using ..SymbolicAnalysis

    function symbolize(e :: Symexpr, tp :: Type)
        if isa(e, Sym)
            if tp <: SparseMatrixCSC
                s = MiddleSymbol(:true)
            else
                dump(tp)
                assert(0)
            end
            return s 
        else
            assert(0)
        end
    end

    function preprocess(property_proxies, symbol_info)
        predefined = Dict{Sym, AbstractSymbol}()
        for (s, v) in property_proxies 
            if isa(v.structure_only, MiddleSymbol)
                assert(v.structure_only.value == :true)
                predefined[s] = v.structure_only
            end
        end 
        predefined
    end

    function postprocess(res, property_proxies, symbol_info)
        for (s, v) in res
            assert(isa(s, Sym))
            if isa(v, MiddleSymbol)
                assert(v.value == :true)
                property_proxies[s].structure_only = v
            end
        end 
    end

    function set_true_action(e)
        e.svalue = MiddleSymbol(:true)
    end

    function pass_a1_action(e)
        e.svalue = e.args[1].svalue
    end

    function assign_action(e)
        e.args[1].svalue = e.args[2].svalue
    end

    const transfer_rules = (
        ((:(=), Any, Any), assign_action),

        ((:call, GlobalRef(Main, :spones), SparseMatrixCSC), set_true_action),
        ((:call, GlobalRef(Main, :speye), Any), set_true_action),

        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Float64), GlobalRef(Main, :Int32)), Any), pass_a1_action),
        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Cdouble), GlobalRef(Main, :Cint)), Any), pass_a1_action),
    )

    const pass_info = ("StructOnly", transfer_rules, preprocess, postprocess, symbolize) 
end
