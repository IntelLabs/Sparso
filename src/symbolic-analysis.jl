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

@doc """
Symbolic data-flow analysis framework 

The framwork is based on 3-level lattice:

    Top
 /  / \  \ 
S1 S2 S3 S4
 \  \ /  /
  Bottom

Each entry in the lattice is a symbol(AbstactSymbol).
"""
module SymbolicAnalysis

using CompilerTools
using SparseAccelerator
using SparseAccelerator: ExprPattern, RegionInfo, AbstractPatternAction, CallSites, expr_skeleton
using SparseAccelerator: dprintln, dprint
import Base.Callable

export TransferRule
export SymbolicAnalyzer
export AbstractSymbol, MiddleSymbol, TOP_SYMBOL, BOTTOM_SYMBOL
export union_symbols

export check_and_gen_transfer_rules
export run_analyzer

export SA_SUCC

typealias BasicBlock      CompilerTools.CFGs.BasicBlock
typealias Statement       CompilerTools.CFGs.TopLevelStatement
typealias Liveness        CompilerTools.LivenessAnalysis.BlockLiveness
typealias CFG             CompilerTools.CFGs.CFG
#typealias DomLoops        CompilerTools.Loops.DomLoops
#typealias Loop            CompilerTools.Loops.Loop
#typealias GenSymId        Int
#typealias BasicBlockIndex Int
#typealias StatementIndex  Int
typealias Sym             Union{Symbol, GenSym} # A Symbol or GenSym.
#typealias Sym2TypeMap     Dict{Sym, Type}
typealias Symexpr         Union{Symbol, GenSym, Expr} # A Symbol, GenSym or Expr

const SA_SUCC = true

abstract AbstractSymbol

type TopSymbol <: AbstractSymbol
end

type BottomSymbol <: AbstractSymbol
end

immutable MiddleSymbol <: AbstractSymbol
    value :: Any
#    class :: SymbolClass

    MiddleSymbol(v) = new(v) #, SymbolClass(v)) 
end

const TOP_SYMBOL = TopSymbol()
const BOTTOM_SYMBOL = BottomSymbol()

@doc """
any ^ Bottom = Bottom
any ^ Top = any 
Mi ^ Mi = Mi
Mi ^ Mj = Bottom (if i != j)
"""
function DefaultMeetRule(s1::AbstractSymbol, s2::AbstractSymbol)
    if s1 == BOTTOM_SYMBOL || s2 == BOTTOM_SYMBOL
        return BOTTOM_SYMBOL
    end

    if s1 == TOP_SYMBOL && s2 != TOP_SYMBOL
        return s2
    elseif s2 == TOP_SYMBOL && s1 != TOP_SYMBOL
        return s1
    elseif s1 == TOP_SYMBOL && s2 == TOP_SYMBOL
        return TOP_SYMBOL
    end

    if s1.value == s2.value #|| s1.class == s2.class
        #assert(isa(s1.class.first, MiddleClass))
        #return s1.class.first
        return s1
    else
        return BOTTOM_SYMBOL
    end
end

function DefaultInitializer()
end

function DefaultSymbolizer(e::Symexpr)
    if isa(e, Sym)
        return MiddleSymbol(Sym)
    else  
        #TODO
    end
end

#function DefaultBooleanSymbolizer(e::Symexpr)
#    return MiddleSymbol(:true)
#end

@doc """
"""
#type SymbolicValue
#    value   :: AbstractSymbol
#    SymbolicValue() = new(TOP_SYMBOL)

#    values  :: Set{AbstractSymbol}  # 
#    final   :: AbstractSymbol       # final symbolic value
 
#    SymbolicValue() = new(Set{AbstractSymbol}(), TOP_SYMBOL)
#end
typealias SymbolicValue AbstractSymbol

@doc """ a default pattern matcher """
function EmptyMatcher(ctx)
    return true
end

@doc """
"""
type TransferRule
    match_pattern   :: Tuple
    match_func      :: Function
    action          :: Function
    
    TransferRule(p :: Tuple, action :: Function) =  new(p, EmptyMatcher, action)
end

type SymbolicAnalyzer
    name            :: AbstractString
    transfer_rules  :: Vector{TransferRule}
    symbolizer      :: Function
    meet_rule       :: Function

    SymbolicAnalyzer(name::AbstractString, 
        transfer_rules::Vector{TransferRule}, 
        symbolizer::Function,
        meet_rule::Function = DefaultMeetRule,
    ) = new(name, transfer_rules, symbolizer, meet_rule)
end

#@doc """
#Context used for analyzing one ast node (Sym or Expr)
#"""
type SymExprContext
    raw_expr    :: Any
    args        :: Vector{SymExprContext}
    types       :: Vector{Type}
    svalue      :: Union{SymbolicValue, Void} #

    SymExprContext(ast) = new(
        ast,
        Vector{SymExprContext}(),
        Vector{Type}(),
        nothing
    )
end


@doc """
Context used for analyzing one Julia statement
"""
type StmtContext
    stmt                 :: Statement                       # statement itself
    in                   :: Set{StmtContext}                # flow-in
    out                  :: Set{StmtContext}                # flow-out

    context_map          :: Dict{Any, SymExprContext}

    changed              :: Set{Sym}                        # variables whose property has been changed?
    # property map after analyzing this statement
    syms                 :: Set{Sym}
    in_symbolic_vals     :: Union{Dict{Sym, SymbolicValue}, Void}        # symbol -> property
    out_symbolic_vals    :: Dict{Sym, SymbolicValue}        # symbol -> property
#    local_symbolic_vals  :: Dict{Symexpr, SymbolicValue}    # symbol|expr -> property

    StmtContext(stmt) = new(
        stmt, # stmt
        Set{StmtContext}(), # in
        Set{StmtContext}(), # out
        #ypeof(stmt.expr) == Expr ? SymExprContext(stmt.expr) : nothing, 
        Dict{Symexpr, SymExprContext}(), # context_map
        Set{Sym}(),                      # changed 
        Set{Sym}(),                      # syms 
        nothing,      # in_symbolic_vals
        Dict{Sym, SymbolicValue}(),      # out_symbolic_vals
#        Dict{Symexpr, SymbolicValue}()   # local_symbolic_vals
    )
end

@doc """
Context used for the whole analysis
"""
type AnalysisContext
    analyzer                :: SymbolicAnalyzer
    stmt_contexts           :: Dict{Statement, StmtContext}
    curr_context            :: Union{StmtContext, Void}

    live_in_before_expr     :: Set{Sym}
    equivalent_sets         :: Dict{SymbolicValue, Vector{SymbolicValue}}
    # Such patterns should be matched at the last, because otherwise, other 
    # patterns may not be able to match what they should: they cannot find the
    # subtrees to match, which are no longer in the same statement.    
    # We call such patterns splitting patterns, and the other non-splitting.
    # We can top-down match non-splitting patterns once, and bottom-up match
    # all patterns (including splitting and non-splitting) once.
    non_splitting_patterns  :: Vector{ExprPattern}

    AnalysisContext(analyzer) = new(
        analyzer,

        Dict{Statement, StmtContext}(), 
        nothing,                         # curr_context

        Set{Sym}(),                      # live_in_before_expr 
        Dict{Sym, Vector{Sym}}(),        # equivalent_sets

        Vector{ExprPattern}()            # non_splitting_patterns
    )
end

function check_and_gen_transfer_rules(rules) 
    ret = Vector{TransferRule}()
    for r in rules
        assert(isa(r, Tuple)) 
        rule = TransferRule(r[1], r[2])
        push!(ret, rule)
    end
    ret
end

type PreActionFunctor <: AbstractPatternAction
    caller :: Function
    action :: Callable
end

function call_pre_action(
    functor           :: PreActionFunctor, 
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
) 
    if typeof(ast) <: Sym || typeof(ast) == SymbolNode 
        ast = SparseAccelerator.get_symexpr(ast)
    end

    stmt_ctx = call_sites.extra.curr_context
    assert(in(ast, stmt_ctx.context_map))
    return functor.action(stmt_ctx.context_map[ast]) # 
end

function generate_pre_action(func::Callable)
    functor = PreActionFunctor(call_pre_action, func)
    functor
end

type PostActionFunctor <: AbstractPatternAction
    caller :: Function
    action :: Callable
end

type PatternActionEnv
    args :: Vector{}
end

function prepare_pattern_action_env(
    stmt_ctx    :: StmtContext,
    expr        :: Expr,
    symbol_info, 
)
    if !in(expr, keys(stmt_ctx.context_map))
        dump(expr)
        dump(keys(stmt_ctx.context_map))
    end
    assert(in(expr, keys(stmt_ctx.context_map)))
    env = stmt_ctx.context_map[expr]
    assert(env.raw_expr == expr)
    
    # fill args if it's empty
    if isempty(env.args)
        idx = env.raw_expr.head == :call ? 2 : 1
        #skeleton = expr_skeleton(expr, symbol_info)
        #dump(skeleton)
        #assert(length(env.raw_expr.args) == length(skeleton))
        while idx <= length(env.raw_expr.args)
            arg = env.raw_expr.args[idx]
            if typeof(arg) == SymbolNode 
                arg = SparseAccelerator.get_symexpr(arg)
            end

            if !in(arg, keys(stmt_ctx.context_map))
                dump(arg)
            end
            assert(in(arg, keys(stmt_ctx.context_map)))
            push!(env.args, stmt_ctx.context_map[arg]) 
            #push!(env.types, skeleton[idx])
            idx = idx + 1
        end
    #else
    #    assert(length(env.args) == length(env.raw_expr.args)-1)
    end

    return env 
end

function done_pattern_action_env(
    stmt_ctx    :: StmtContext,
    env         :: SymExprContext
)
    #stmt_ctx.syms
    #stmt_ctx.out_symbolic_vals
end


function call_post_action(
    functor           :: PostActionFunctor, 
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
) 
    if typeof(ast) <: Sym || typeof(ast) == SymbolNode 
        ast = SparseAccelerator.get_symexpr(ast)
    end

    # wrap env
    stmt_ctx = call_sites.extra.curr_context
    env = prepare_pattern_action_env(stmt_ctx, ast, call_sites.symbol_info)
    ret = functor.action(env) # 
    # assign a unique symbolic value if expr matched but not processed 
    # by the user defined pattern rules
    #if ret != SA_SUCC && env.svalue == nothing && isa(env.raw_expr, Symexpr) && env.raw_expr.head != :(=) 
    #    env.svalue = call_sites.extra.analyzer.symbolizer(env.raw_expr, Any, true)
    #end
    done_pattern_action_env(stmt_ctx, env)
    return true
end

function generate_post_action(func::Callable)
    functor = PostActionFunctor(call_post_action, func)
    functor
end

function generate_analysis_patterns(
    analyzer :: SymbolicAnalyzer
)
    # generate propagation patterns from user defined rules (analyzer)
    analysis_patterns = Vector{ExprPattern}()
    for r in analyzer.transfer_rules
        assert(isa(r.match_pattern, Tuple))
        pre_action = r.match_func == EmptyMatcher ? SparseAccelerator.do_nothing : generate_pre_action(r.match_func)
        post_action = generate_post_action(r.action) 
        #post_action =  SparseAccelerator.do_nothing
        push!(analysis_patterns, 
            ExprPattern(string(r), # name
                r.match_pattern,
                (:NO_SUB_PATTERNS,),
                pre_action,
                (:NO_CHANGE,),
                post_action,  
                "",
                "",
                (),
                0, 
                ()
            )
        )
    end
    return analysis_patterns
end

@doc """ 
build context for every top-level statement in a code region
"""
function build_stmt_contexts(
    region :: RegionInfo
) 
    stmt_to_ctx = Dict{Statement, StmtContext}()
    for s in region.stmts
        #dprintln(1, 0, s)
        stmt_to_ctx[s] = StmtContext(s)
    end

    bb_idxs = keys(region.bblocks)
    bbs = values(region.bblocks)
    for bb in bbs
        if isempty(bb.statements)
            continue
        end

        first_stmt = first(bb.statements)
        last_stmt = last(bb.statements)

        # intrablock in and out
        for (idx, s) in enumerate(bb.statements)
            if s == last_stmt # last one?
                break
            end
            next = bb.statements[idx+1]
            assert(in(s, keys(stmt_to_ctx)))
            assert(in(next, keys(stmt_to_ctx)))
            push!(stmt_to_ctx[s].out, stmt_to_ctx[next])
            push!(stmt_to_ctx[next].in, stmt_to_ctx[s])
        end

        # interblock in and out
        # backward
        bb_list = collect(copy(bb.preds))
        i = 1
        while i <= length(bb_list)
            pred = bb_list[i]
            i+=1
            if !in(pred, bbs)
                continue
            end
            if isempty(pred.statements)
                for p in pred.preds
                    if !in(p, bb_list)
                        push!(bb_list, p)
                    end
                end
            else
                pred_last = last(pred.statements)
                push!(stmt_to_ctx[first_stmt].in, stmt_to_ctx[pred_last])
            end
        end
        # forward
        bb_list = collect(copy(bb.succs))
        i = 1
        while i <= length(bb_list)
            succ = bb_list[i]
            i+=1
            if !in(succ, bbs)
                continue
            end
            if isempty(succ.statements)
                for p in succ.succs
                    if !in(p, bb_list)
                        push!(bb_list, p)
                    end
                end
            else
                succ_first = first(succ.statements)
                push!(stmt_to_ctx[last_stmt].out, stmt_to_ctx[succ_first])
            end
        end
    end

    return stmt_to_ctx
end

# Symbol types unimportant to analysis
const SKIP_TYPES = [GlobalRef, TopNode, Int32, Int64, Float64, Bool, QuoteNode, ASCIIString, Complex]

function build_expr_contexts(ast, stmt_ctx)
    expr = ast
    ast_type = typeof(ast)
    if ast_type <: Sym || ast_type <: SymbolNode 
        expr = SparseAccelerator.get_symexpr(ast)        
        ctx = SymExprContext(expr)
        ctx.svalue = TOP_SYMBOL
        stmt_ctx.context_map[expr] = ctx
        push!(stmt_ctx.syms, expr)
    elseif ast_type <: Expr
        ctx = SymExprContext(expr)
        #ctx.svalue = TOP_SYMBOL
        stmt_ctx.context_map[expr] = ctx
    #elseif isa(expr, TopNode)
        #dump(expr)
        #assert(0)
    #elseif isa(expr, GlobalRef)
    elseif in(ast_type, SKIP_TYPES)
        #expr = SparseAccelerator.get_symexpr(ast)        
        #dump(ast)
        ctx = SymExprContext(expr)
        #ctx.svalue = TOP_SYMBOL
        stmt_ctx.context_map[expr] = ctx
    else
        assert(0)
        #return nothing
    end

    #if !in(expr, keys(stmt_ctx.context_map))
    #    ctx = SymExprContext(expr)
    #    stmt_ctx.context_map[expr] = ctx
    #    if isa(ast, Expr)
    #        for i in 2 : length(expr.args)
    #            push!(ctx.args, (build_expr_contexts(ast.args[i], stmt_ctx)))
    #       end
    #    end
    #    return ctx
    #end
    return nothing
end


function build_expr_contexts(ast, stmt_ctx, top_level_number, is_top_level, read)
    build_expr_contexts(ast, stmt_ctx) 
end

function merge_symbolic_vals(des::Dict{Sym, SymbolicValue}, src::Dict{Sym, SymbolicValue}, meet_rule::Function)
    des_keys = collect(keys(des))
    for (k, v) in src
        if in(k, des_keys)
            des[k] = meet_rule(src[k], des[k])
        else
            des[k] = v
        end
    end
end

#function show_
#end

function show_contexts(stmt_contexts::Dict{Statement, StmtContext})
    for (s, ctx) in stmt_contexts
        dprintln(1, 0, s.index, " ", s.expr)
        dprint(1, 1, "IN:")
        for c in ctx.in
            dprint(1, 0, " ", c.stmt.index)
        end
        dprintln(1, 0, "")
        dprint(1, 1, "OUT:")
        for c in ctx.out
            dprint(1, 0, " ", c.stmt.index)
            #dprintln(1, 2, c.stmt.index, " ", c.stmt.expr)
        end
        dprintln(1, 0, "")
    end
end

function union_symbols(s1::MiddleSymbol, s2::MiddleSymbol)
    if s1.value != s2.value && s1.class != s2.class
        assert(isempty(intersect(s1.class, s2.class)))
        union!(s1.class, s2.class)
        s2.class = s1
    end
end

function show_symbolic_vals(vals::Dict{Sym, SymbolicValue}, indent_level=1)
    for (s, v) in vals
        if v == BOTTOM_SYMBOL
            val = "BOTTOM"
        elseif v == TOP_SYMBOL
            val = "TOP"
        else
            val = v.value
        end

        dprintln(1, indent_level, s, " -> ", val)
        # meeting func
    end
end

function run_analyzer(
    analyzer    :: SymbolicAnalyzer,
    region      :: RegionInfo,
    predefined  ::Dict{Sym, AbstractSymbol}, 
    verbose 
) 

    # initialize property values for all Syms according to

    analysis_context = AnalysisContext(analyzer)

    analysis_context.stmt_contexts = build_stmt_contexts(region) 

    # show_contexts(analysis_context.stmt_contexts)

    # build contexts for expressions of every statement
    for (stmt, ctx) in analysis_context.stmt_contexts
        expr = stmt.expr
        if typeof(expr) != Expr || expr.head == :return
            continue
        end
        CompilerTools.AstWalker.AstWalk(expr, build_expr_contexts, ctx)
    end

    analysis_patterns = generate_analysis_patterns(analyzer)

    call_sites  = SparseAccelerator.CallSites(Set{SparseAccelerator.CallSite}(), 
                            region.region, 
                            region.lambda, 
                            region.symbol_info,
                            region.liveness, analysis_patterns,
                            Vector{SparseAccelerator.Action}(), analysis_context)

  
    # All patterns are non-splittable.
    analysis_context.non_splitting_patterns = analysis_patterns
 
    # Build flow graph

    converged = false
    cnt = 0
    while !converged
        cnt = cnt + 1
        converged = true
        for (stmt, stmt_ctx) in analysis_context.stmt_contexts

            if verbose > 1
                dprintln(1, 0, stmt.index, " ", stmt.expr)
            end

            new_symbolic_vals = Dict{Sym, SymbolicValue}()
            for in_ctx in stmt_ctx.in
                merge_symbolic_vals(new_symbolic_vals, in_ctx.out_symbolic_vals, analyzer.meet_rule)
            end

            # preserve predefined values
            for s in keys(predefined)
                #if in(s, keys(stmt_ctx.in_symbolic_vals))
                new_symbolic_vals[s] = predefined[s] #analyzer.symbolizer(s)
                #end
            end

            if verbose > 1
                dprintln(1, 1, "IN ")
                show_symbolic_vals(new_symbolic_vals, 2)
            end

            expr = stmt.expr
            if stmt_ctx.in_symbolic_vals == new_symbolic_vals 
                if verbose > 1
                    dprintln(1, 1, "IN_SET_NOT_CHANGED")
                end
                continue
            else
                stmt_ctx.in_symbolic_vals = new_symbolic_vals
            end

            if typeof(expr) != Expr || expr.head == :return
                stmt_ctx.out_symbolic_vals = new_symbolic_vals
                if verbose > 1
                    dprintln(1, 1, "SKIP")
                    show_symbolic_vals(new_symbolic_vals, 2)
                end
                continue
            end

            #ctx = analysis_context.stmt_contexts[stmt]
            #ctx.local_map.clear()
            analysis_context.curr_context = stmt_ctx
            call_sites.extra.live_in_before_expr = CompilerTools.LivenessAnalysis.live_in(stmt, region.liveness)

            # update all Symexpr contexts from new_symbolic_vals
            for (se, se_ctx) in stmt_ctx.context_map
                if in(se, keys(new_symbolic_vals)) #(isa(se, Symbol) || isa(se, SymbolNode))
                    se_ctx.svalue = new_symbolic_vals[se]
                end
            end

            CompilerTools.AstWalker.AstWalk(expr, SparseAccelerator.match_replace_an_expr_pattern, call_sites)
 
            # update out_symbolic_vals
            for (se, se_ctx) in stmt_ctx.context_map
                if isa(se, Sym) #&& se_ctx.svalue != nothing 
                    if se_ctx.svalue == nothing
                        #dump(se_ctx)
                        #assert(0)
                        se_ctx.svalue = TOP_SYMBOL
                    end
                    assert(in(se, stmt_ctx.syms))
                    if !in(se, keys(new_symbolic_vals))
                        new_symbolic_vals[se] = se_ctx.svalue
                    elseif se_ctx.svalue != new_symbolic_vals[se]
                        #if se_ctx.svalue == nothing
                        #end
                        new_symbolic_vals[se] = se_ctx.svalue
                    end
                end
            end

            # preserve predefined values
            for s in keys(predefined)
                #if in(s, keys(stmt_ctx.out_symbolic_vals))
                new_symbolic_vals[s] = predefined[s] #analyzer.symbolizer(s)
                #end
            end

            if stmt_ctx.out_symbolic_vals != new_symbolic_vals
                stmt_ctx.out_symbolic_vals = new_symbolic_vals
                converged = false

                if verbose > 1
                    dprintln(1, 1, "OUT")
                    show_symbolic_vals(new_symbolic_vals, 2)
                end
            else
                if verbose > 1
                    dprintln(1, 1, "OUT_SET_NOT_CHANGED")
                    show_symbolic_vals(new_symbolic_vals, 2)
                end
            end
        end
    end 

    # global symbolic values
    global_symbolic_vals = Dict{Sym, SymbolicValue}()
    for (stmt, stmt_ctx) in analysis_context.stmt_contexts
        #merge_symbolic_vals(global_symbolic_vals, stmt_ctx.in_symbolic_vals, analyzer.meet_rule)
        merge_symbolic_vals(global_symbolic_vals, stmt_ctx.out_symbolic_vals, analyzer.meet_rule)
    end

    if verbose > 0
        dprintln(1, 0, "\n", region.name, " ", analyzer.name)
        show_symbolic_vals(global_symbolic_vals)
    end

    return global_symbolic_vals
end

end # module
