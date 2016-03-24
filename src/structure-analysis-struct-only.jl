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

module StructureAnalysisStructOnly
   
    using SparseAccelerator: Sym, Symexpr, TypedExprNode
    using ..SymbolicAnalysis

    function symbolize(e :: Symexpr, tp :: Type, unique)
        if unique == true 
            sym = new_symbol("unknown")
            return MiddleSymbol(sym)
        end 

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
            if isa(v, MiddleSymbol) && v.value == :true
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

    function set_elem_action(e)
        dump(e.args[2].raw_expr)
        val = e.args[2].raw_expr
        if !isa(val, Int) || (val != 0 && val != 1)
            e.args[1].svalue = MiddleSymbol(:false)
        end
    end

    const transfer_rules = (
        ((:(=), Any, Any), assign_action),

        ((:call, GlobalRef(Main, :spones), Any), set_true_action),
        ((:call, GlobalRef(Main, :zeros), Any, Any), set_true_action),
        ((:call, GlobalRef(Main, :sparse), Any), pass_a1_action),
        ((:call, GlobalRef(Main, :scale), SparseMatrixCSC, Any), pass_a1_action),
        ((:call, GlobalRef(Main, :setindex!), Array, Any, Any, Any), set_elem_action),

        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Float64), GlobalRef(Main, :Int32)), Any), pass_a1_action),
        ((:call, TypedExprNode(Function, :call, TopNode(:apply_type), GlobalRef(Main, :SparseMatrixCSC), GlobalRef(Main, :Cdouble), GlobalRef(Main, :Cint)), Any), pass_a1_action),
    )

    const pass_info = ("StructOnly", transfer_rules, preprocess, postprocess, symbolize) 
end
