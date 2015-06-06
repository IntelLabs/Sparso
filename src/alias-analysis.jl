# A basic alias analysis with limited scope:
#
# 1. Only objects pointed to by variables, but not elements of arrays or other structs
# 2. We only consider variables that are assigned once
# 3. No inter-precedural (i.e. function call) analysis
#
# The only useful result from this alias analysis is whether
# some variable definitely doesn't alias with anything else.
#
# We are NOT interested in the set of potential aliased variables.
#
# The algorithm is basically an abstract interpreter of Julia AST.

module AliasAnalysis

using CompilerTools

DEBUG_LVL=0

function set_debug_level(x)
    global DEBUG_LVL = x
end

# A debug print routine.
function dprint(level,msgs...)
    if(DEBUG_LVL >= level)
        print(msgs...)
    end 
end

# A debug print routine.
function dprintln(level,msgs...)
    if(DEBUG_LVL >= level)
        println(msgs...)
    end 
end

# state to keep track of variable values
const Unknown  = -1
const NotArray = 0

type State
  baseID :: Int
  locals :: Dict{Symbol, Int}
  revmap :: Dict{Int, Set{Symbol}}
  nest_level :: Int
  top_level_idx :: Int
  liveness :: CompilerTools.LivenessAnalysis.BlockLiveness
end

init_state(liveness) = State(0, Dict{Symbol,Int}(), Dict{Int, Set{Symbol}}(), 0, 0, liveness)

function next_node(state)
  local n = state.baseID + 1
  state.baseID = n
  return n
end

function update_node(state, v, w)
  if !haskey(state.locals, v)
    # new initialization
    push!(state.locals, v, w)
    if haskey(state.revmap, w)
      push!(state.revmap, w, push!(state.revmap[w], v))
    else
      push!(state.revmap, w, push!(Set{Symbol}(), v))
    end
  else
    # when a variable is initialized more than once, set to Unknown
    push!(state.locals, v, Unknown)
    if haskey(state.revmap, w)
      for u in state.revmap[w]
        push!(state.locals, u, Unknown)
      end
      pop!(state.revmap, w)
    end
  end
end

function update_unknown(state, v)
  state.locals[v] = Unknown
end

function update_notarray(state, v)
  state.locals[v] = NotArray
end

function update(state, v, w)
  if w > 0
    update_node(state, v, w)
  elseif w == NotArray
    update_notarray(state, v)
  else
    update_unknown(state, v)
  end
end

function lookup(state, v)
  if haskey(state.locals, v)
    state.locals[v]
  else
    Unknown
  end
end

# (:lambda, {param, meta@{localvars, types, freevars}, body})
function from_lambda(state, env, expr)
  local head = expr.head
  local ast  = expr.args
  local typ  = expr.typ
  assert(length(ast) == 3)
  local param = ast[1]
  local meta  = ast[2] # { {Symbol}, {{Symbol,Type,Int}}, {Symbol,Type,Int} }
  local body  = ast[3]
  # very conservative handling by setting free variables to Unknown.
  # TODO: may want to simulate function call at call site to get 
  #       more accurate information.
  for (v,vt,vm) in meta[3]
    update_unknown(state, v)
  end
  return NotArray
end

function from_exprs(state, env, ast)
  local len  = length(ast)
  [ from_expr(state, env, exp) for exp in ast ]
end
 
function from_body(state, env, expr::Any)
  local exprs = expr.args
  local ret = NotArray       # default return
  for i = 1:length(exprs)
    if state.nest_level == 0
      state.top_level_idx = i
    end
    ret = from_expr(state, env, exprs[i])
  end
  return ret
end

function from_assignment(state, env, expr::Any) 
  local head = expr.head
  local ast  = expr.args
  local typ  = expr.typ
  assert(length(ast) == 2)
  local lhs = ast[1]
  local rhs = ast[2]
  dprintln(2, "AA ", lhs, " = ", rhs)
  if (isa(lhs, SymbolNode)) 
    lhs = lhs.name
  end
  assert(isa(lhs, Symbol))
  if lookup(state, lhs) != NotArray
    rhs = from_expr(state, env, rhs)
    # if all vars that have rhs are not alive afterwards
    # then we can safely give v a fresh ID.
    if state.nest_level == 0 
      tls = LivenessAnalysis.find_top_number(state.top_level_idx, state.liveness)
      assert(tls != nothing)
      assert(LivenessAnalysis.isDef(lhs, tls))
      if (haskey(state.revmap, rhs))
        dead = true
        for v in state.revmap[rhs]
          dead = dead && !in(v, tls.live_out)
        end
        if dead
          rhs = next_node(state)
        end
      end
    end
    dprintln(2, "AA update ", lhs, " <- ", rhs)
    update(state, lhs, rhs)
  end
end

function normalize_callname(state, env, fun, args)
  if isa(fun, Symbol)
    if is(fun, :broadcast!)
      dst = args[2]
      if isa(dst, SymbolNode)
        dst = get(state.defs, dst.name, nothing)
        if isa(dst, VarDef)
          dst = dst.rhs
        end
      end
      if isa(dst, Expr) && is(dst.head, :call) && isa(dst.args[1], TopNode) &&
         is(dst.args[1].name, :ccall) && isa(dst.args[2], QuoteNode) &&
         is(dst.args[2].value, :jl_new_array)
        # now we are sure destination array is temporary
        fun   = args[1]
        args  = args[3:end]
        if isa(fun, Symbol)
        elseif isa(fun, SymbolNode)
          fun = fun.name
        else
          error("cannot handle broadcast! with function ", fun)
        end
      else
        dprintln(env, "cannot decide :broadcast! destination is temporary")
      end
    end
  elseif isa(fun, TopNode)
    fun = fun.name
    if is(fun, :ccall)
      callee = args[1]
      if isa(callee, SymbolNode)
        callee = get(state.defs, callee.name, nothing)
      end
      if isa(callee, QuoteNode) && (is(callee.value, :jl_alloc_array_1d) || is(callee.value, :jl_alloc_array_2d) || is(callee.value, :jl_alloc_array_3d))
        local realArgs = Any[]
        for i = 4:2:length(args)
          push!(realArgs, args[i])
        end
        fun  = :alloc
        args = realArgs
      else
      end
    end
  elseif isa(fun, GlobalRef)
    if is(fun.mod, Base.Broadcast)
      if is(fun.name, :broadcast_shape)
        fun = :broadcast_shape
      end
    end
  end
  return (fun, args)
end

mapSym = Symbol[:negate, :.+, :.-, :.*, :./, :+, :-, :*, :/, :sin, :erf, :log10, :exp, :sqrt]

function from_call(state, env, expr::Any)
  local head = expr.head
  local ast = expr.args
  local typ = expr.typ
  assert(length(ast) >= 1)
  local fun  = ast[1]
  local args = ast[2:end]
  dprintln(2, "AA from_call: fun=", fun, " typeof(fun)=", typeof(fun), " args=",args, " typ=", typ)
  #fun = from_expr(state, env, fun)
  #dprintln(2, "AA from_call: new fun=", fun)
  (fun_, args) = normalize_callname(state, env, fun, args)
  dprintln(2, "AA from_call: normalized fun=", fun_)
  if is(fun_, :arrayref) || is(fun, :arrayset) || in(fun_, mapSym)
    # This is actually an conservative answer since arrayref might return
    # an array too, but we don't consider it as a case to handle.
    return NotArray
  elseif is(fun_, :alloc)
    return next_node(state)
  elseif is(fun_, :fill!)
    return from_expr(state, env, args[1])
  else
    dprintln(2, "AA: unknown call ", fun_)
    # For unknown calls, conservative assumption is that after 
    # the call, its array type arguments might alias each other.
    for exp in args
      if isa(exp, SymbolNode)
        update_unknown(state, exp.name)
      elseif isa(exp, Symbol)
        update_unknown(state, exp)
      end
    end
    return Unknown
  end
end

function from_return(state, env, expr)
  local head = expr.head
  local typ  = expr.typ
  local args = from_exprs(state, env, expr.args)
  if length(args) == 1
    return args[1]
  else
    return Unknown
  end
end

function from_expr(state, env, ast)
  if isa(ast, LambdaStaticData)
      ast = uncompressed_ast(ast)
  end
  local asttyp = typeof(ast)
  dprint(2, "AA from_expr: ", asttyp)
  if is(asttyp, Expr)
    local head = ast.head
    local args = ast.args
    local typ  = ast.typ
    dprintln(2, " --> ", head)
    if is(head, :lambda)
        return from_lambda(state, env, ast)
    elseif is(head, :body)
        return from_body(state, env, ast)
    elseif is(head, :(=))
        return from_assignment(state, env, ast)
    elseif is(head, :return)
        return from_return(state, env, ast)
    elseif is(head, :call)
        return from_call(state, env, ast)
        # TODO: catch domain IR result here
    elseif is(head, :call1)
        return from_call(state, env, ast)
    elseif is(head, :assertEqShape)
        return NotArray
    elseif is(head, :alloc)
        return next_node(state)
    elseif is(head, :copy)
        return next_node(state)
    elseif is(head, :method)
        # skip
    elseif is(head, :line)
        # skip
    elseif is(head, :new)
        # skip
    elseif is(head, :gotoifnot)
        # skip
    elseif is(head, :loophead)
        # skip
    elseif is(head, :loopend)
        # skip
    else
        throw(string("from_expr: unknown Expr head :", head))
    end
  elseif is(asttyp, SymbolNode)
    dprintln(2, " ", ast)
    return lookup(state, ast.name)
  else
    dprintln(2, " not handled ", ast)
  end
  return Unknown
end

function isarray(typ)
  isa(typ, DataType) && is(typ.name, Array.name)
end

# (:lambda, {param, meta@{localvars, types, freevars}, body})
function analyze_lambda_body(body, param, meta_typed, liveness)
  local state = init_state(liveness)
  dprintln(2, "AA ", isa(body, Expr), " ", is(body.head, :body)) 
  # FIXME: surprisingly the first value printed above is false!
  for (v, (u,t,m)) in meta_typed
    if !isarray(t)
      update_notarray(state, v)
    end
  end
  for v in param
    # Note we assume all input parameters do not aliasing each other,
    # which is a very strong assumption. This may require reconsideration.
    # Update: changed to assum nothing by default.
    if isarray(meta_typed[v][2])
      #update_node(state, v, next_node(state))
      update_unknown(state, v)
    end
  end
  dprintln(2, "AA locals=", state.locals)
  from_expr(state, Nothing, body)
  dprintln(2, "AA locals=", state.locals)
  local revmap = Dict{Int, Symbol}()
  local unique = Set{Symbol}()
  # keep only variables that have unique object IDs. 
  # TODO: should consider liveness either here or during analysis, 
  #       since its ok to alias dead vars.
  for (v, w) in state.locals
    if w > 0
      if haskey(revmap, w)
        delete!(unique, revmap[w])
      else
        push!(unique, v)
        push!(revmap, w, v)
      end
    end
  end
  dprintln(2, "AA after alias analysis: ", unique)
  # return the set of variables that are confirmed to have no aliasing
  return unique  
end

# (:lambda, {param, meta@{localvars, types, freevars}, body})
function analyze_lambda(expr, liveness)
  local head = expr.head
  local ast  = expr.args
  assert(length(ast) == 3)
  local param = ast[1]
  local meta  = ast[2] # { {Symbol}, {{Symbol,Type,Int}}, {Symbol,Type,Int} }
  local body  = ast[3]
  local meta2 = Dict{Symbol, Array{Any}}()
  for info in meta[2]
    meta2[info[1]] = info
  end
#  return nothing
  analyze_lambda_body(body, param, meta2, liveness)
end

end

