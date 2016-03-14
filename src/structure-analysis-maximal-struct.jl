#=
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

module StructureAnalysisMaximalStruct
    
    import SparseAccelerator
    using SparseAccelerator: Sym, Symexpr, TypedExprNode
    using ..SymbolicAnalysis
    using ..SymbolicAnalysis: SA_SUCC

    function symbolize(e :: Symexpr, tp :: Type, unique = false)
        if isa(e, Sym)
            if tp <: Vector
                s = MiddleSymbol((e))
            elseif tp <: SparseMatrixCSC
                s = MiddleSymbol((e))
            elseif tp <:Int || tp <: Float64
                s = MiddleSymbol((:NUM_1))
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
            if !isa(s, Symbol)
                continue
            end
            if isa(v.constant_valued, MiddleSymbol)
                predefined[s] = symbolize(s, symbol_info[s]) 
            end

            if isa(v.maximal_structured, MiddleSymbol)
                predefined[s] = v.maximal_structured
            end
        end 
        predefined
    end

    function postprocess(res, property_proxies, symbol_info)
        for (s, v) in res
            assert(isa(s, Sym))
            if isa(v, MiddleSymbol)
                property_proxies[s].maximal_structured = v
            end
        end 
    end

    function add_sub_action(e)
        if all(a -> a.svalue == BOTTOM_SYMBOL, e.args)
            e.svalue = BOTTOM_SYMBOL
        else
            for arg in e.args
                if isa(arg.svalue, MiddleSymbol)
                    e.svalue = arg.svalue
                    break
                end
            end
        end
    end

    # A * B -> (A, B)
    function A_mul_B_action(e)
        if isa(e.args[1].svalue, MiddleSymbol) && isa(e.args[2].svalue, MiddleSymbol) 
            e.svalue = MiddleSymbol((e.args[1].svalue.value[1], e.args[2].svalue.value[2]))
        end
    end

    # A * B * C
    function mul3_action(e)
        if isa(e.args[1].svalue, MiddleSymbol) && isa(e.args[3].svalue, MiddleSymbol) 
            e.svalue = MiddleSymbol((e.args[1].svalue.value[1], e.args[3].svalue.value[2]))
            if !isa(e.args[2].svalue, MiddleSymbol)
                e.args[2].svalue = e.svalue
            end
        end
    end

    ## v * M -> size of v
    function v_mul_M_action(e)
        assert(0)
        e.svalue = e.args[1].svalue
    end

    # M * v -> (M, M) or (M, v) ?
    function M_mul_v_action(e)
        #dump(e)
        #assert(0)
        if isa(e.args[1].svalue, MiddleSymbol)
            e.svalue = MiddleSymbol((e.args[1].svalue.value[1], :NUM_1))
        end
    end

    function elem_mul_action(e)
        if isa(e.args[1].svalue, MiddleSymbol)
            e.svalue = e.args[1].svalue
        elseif isa(e.args[2].svalue, MiddleSymbol)
            e.svalue = e.args[2].svalue
        end
    end

    function speye_action(e)
        e.svalue = MiddleSymbol((Symbol(e.args[1].raw_expr), Symbol(e.args[1].raw_expr)))
        dump(e.svalue)
    end

    function transpose_action(e)
        if isa(e.args[1].svalue, MiddleSymbol)
            e.svalue = MiddleSymbol((e.args[1].svalue.value[2], e.args[1].svalue.value[1]))
        end
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
        ((:(=), Any, Vector), assign_action),

        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC), A_mul_B_action),
        ((:call, GlobalRef(Main, :spmatmul_witheps), SparseMatrixCSC, SparseMatrixCSC, Any), A_mul_B_action),
        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC), mul3_action),
        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector), M_mul_v_action),
        ((:call, GlobalRef(Main, :*), Vector, SparseMatrixCSC), v_mul_M_action),
        ((:call, GlobalRef(Main, :+), SparseMatrixCSC, SparseMatrixCSC), add_sub_action),
        ((:call, GlobalRef(Main, :+), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC), add_sub_action),
        ((:call, GlobalRef(Main, :-), Any, SparseMatrixCSC), add_sub_action),
        ((:call, GlobalRef(Main, :ctranspose), SparseMatrixCSC), transpose_action),
        ((:call, GlobalRef(Main, :cholfact_int32), SparseMatrixCSC), pass_a1_action),
        ((:call, GlobalRef(Main, :speye), Int), speye_action),

        ((:call, GlobalRef(Main, :*), SparseMatrixCSC, Union{Float64, Int64, Int32}), pass_a1_action),
        ((:call, GlobalRef(Main, :/), SparseMatrixCSC, Union{Float64, Int64, Int32}), pass_a1_action),
        ((:call, GlobalRef(Main, :*), Union{Float64, Int64, Int32}, SparseMatrixCSC), pass_a2_action),

        ((:call, GlobalRef(Main, :.*), Vector, Vector), elem_mul_action),
        ((:call, GlobalRef(Main, :+), Vector, Vector), add_sub_action),
        ((:call, GlobalRef(Main, :-), Vector, Vector), add_sub_action),
        ((:call, GlobalRef(Main, :-), Vector, Vector), add_sub_action),

        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Float64), GlobalRef(Main, :Int32)), Any), pass_a1_action),
    )

    const empty_rules = ()

    const pass_info = ("MaximalStruct", empty_rules, preprocess, postprocess, symbolize) 
end
