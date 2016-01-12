module StructureAnalysisSymmValue
    
    using SparseAccelerator: Sym, Symexpr, TypedExprNode
    using ..SymbolicAnalysis

    function symbolize(e :: Symexpr, tp :: Type)
        if isa(e, Sym)
            if tp <: Vector
                s = MiddleSymbol(:true)
            elseif tp <: SparseMatrixCSC
                s = MiddleSymbol(:true)
            elseif tp <:Int || tp <: Float64
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
        pproxies = property_proxies
        predefined = Dict{Sym, AbstractSymbol}()
        for (s, v) in property_proxies 
            if !isa(s, Symbol)
                continue
            end
            if isa(v.symmetric_valued, MiddleSymbol)
                assert(v.symmetric_valued.value == :true)
                predefined[s] = v.symmetric_valued
            end
        end 
        predefined
    end

    function postprocess(res, property_proxies, symbol_info)
        for (s, v) in res
            assert(isa(s, Symbol))
            if isa(v, MiddleSymbol)
                assert(v.value == :true)
                property_proxies[s].symmetric_valued = v
            end
        end 
    end

    function add_sub_action(e)
        if any(a -> a.svalue == BOTTOM_SYMBOL, e.args)
            e.svalue = BOTTOM_SYMBOL
        elseif all(a -> isa(a.svalue, MiddleSymbol), e.args)
            e.svalue = MiddleSymbol(:true)
        end
    end

    # A * B
    function A_mul_B_action(e)
        A = e.args[1].raw_expr
        B = e.args[2].raw_expr
        if isa(A, SymbolNode) && typeof(B) <: Expr && B.head == :call
            m, func_name = resolve_call_names(B.args)
            if func_name == "ctranspose" && B.args[2].name == A.name
                e.svalue = MiddleSymbol(:true)
            end
        end
    end

    function A_mul_Bc_action(e)
        A = e.args[1].raw_expr
        B = e.args[2].raw_expr
        if isa(A, SymbolNode) && isa(B, SymbolNode) && A.name == B.name
            e.svalue = MiddleSymbol(:true)
        end
    end

    function speye_action(e)
        e.svalue = MiddleSymbol(:true)
    end

    function pass_a1_action(e)
        #dump(e.args[1])
        e.svalue = e.args[1].svalue
    end

    function pass_a2_action(e)
        e.svalue = e.args[2].svalue
    end

    function assign_action(e)
        e.args[1].svalue = e.args[2].svalue
    end

    const transfer_rules = (
        ((:(=), Any, SparseMatrixCSC), assign_action),
        ((:(=), Factorization, Any), assign_action),

        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC), A_mul_B_action),
        ((:call, GlobalRef(Main, :spmatmul_witheps), SparseMatrixCSC, SparseMatrixCSC, Any), A_mul_B_action),

        ((:call, GlobalRef(Main, :A_mul_Bc), SparseMatrixCSC, SparseMatrixCSC), A_mul_Bc_action),
#        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC), mul3_action),
#        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector), M_mul_v_action),
#        ((:call, GlobalRef(Main, :*), Vector, SparseMatrixCSC), v_mul_M_action),
        ((:call, GlobalRef(Main, :+), SparseMatrixCSC, SparseMatrixCSC), add_sub_action),
        ((:call, GlobalRef(Main, :+), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC), add_sub_action),
        ((:call, GlobalRef(Main, :-), Any, SparseMatrixCSC), add_sub_action),
#        ((:call, GlobalRef(Main, :cholfact_int32), SparseMatrixCSC), pass_a1_action),
        ((:call, GlobalRef(Main, :speye), Int), speye_action),

        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, Union{Float64, Int64, Int32}), pass_a1_action),
        ((:call, GlobalRef(Main, :/), SparseMatrixCSC, Union{Float64, Int64, Int32}), pass_a1_action),
        ((:call, GlobalRef(Main, :*), Union{Float64, Int64, Int32}, SparseMatrixCSC), pass_a2_action),

        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Float64), GlobalRef(Main, :Int32)), Any), pass_a1_action),
    )

    const pass_info = ("SymmetricValue", transfer_rules, preprocess, postprocess, symbolize) 
end
