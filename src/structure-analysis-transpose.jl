module StructureAnalysisTranspose
    
    using SparseAccelerator: Sym, Symexpr, TypedExprNode
    using ..SymbolicAnalysis

    function symbolize(e :: Symexpr, tp :: Type, unique = false)
        if unique == true 
            sym = new_symbol("unknown")
            return MiddleSymbol(sym)
        end 

        if isa(e, Sym)
            if tp <: SparseMatrixCSC
                s = MiddleSymbol(e)
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
            if isa(v.transpose_of, MiddleSymbol)
                assert(isa(v.transpose_of.value, Sym))
                predefined[s] = v.transpose_of
            end
        end 
        predefined
    end

    function postprocess(res, property_proxies, symbol_info)
        for (s, v) in res
            assert(isa(s, Sym))
            if isa(v, MiddleSymbol)
                assert(isa(v.value, Sym))
                property_proxies[s].transpose_of = v
            end
        end 
    end

    function transpose_action(e)
        e.svalue = symbolize(e.args[1].raw_expr, SparseMatrixCSC)
    end

    function pass_a1_action(e)
        #dump(e.args[1])
        e.svalue = e.args[1].svalue
    end

    function assign_action(e)
        e.args[1].svalue = e.args[2].svalue
    end

    const transfer_rules = (
        ((:(=), Any, SparseMatrixCSC), assign_action),
        ((:(=), Factorization, Any), assign_action),
        ((:(=), Any, Vector), assign_action),

        ((:call, GlobalRef(Main, :ctranspose), SparseMatrixCSC), transpose_action),

        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Float64), GlobalRef(Main, :Int32)), Any), pass_a1_action),
    )

    const pass_info = ("Transpose", transfer_rules, preprocess, postprocess, symbolize) 
end
