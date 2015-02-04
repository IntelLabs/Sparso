module SparseAccelerator

export @acc

include("ast_walk.jl")
include("liveness.jl")
include("alias-analysis.jl")

#using ..AstWalker
#using ..LivenessAnalysis
#using ..AliasAnalysis

# This controls the debug print level.  0 prints nothing.  At the moment, 2 prints everything.
DEBUG_LVL=3

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

include("sparse-analyze.jl")

function typeOfOpr(x)
  if isa(x, Expr) x.typ
  elseif isa(x, SymbolNode) x.typ
  else typeof(x) 
  end
end   

type memoizeState
  mapNameFuncInfo :: Dict{Any, Any}

  function memoizeState()
    new (Dict{Any,Any}())
  end
end

function findInvariants(members::Set, uniqSet::Set{Symbol}, bbs)
    invariants = Set{Symbol}()

    # Add all the variables that are read and are unique.
    for bbnum in members
      bb = bbs[bbnum]
      for use_var in bb.use
        if in(uniqueSet, use_var)
          push!(invariants, use_var)
        end        
      end 
    end

    # Remove all the variables that may be read in one basic block but written in another.
    for bbnum in members
      bb = bbs[bbnum]

      for def_var in bb.def
        delete!(invariants, def_var) 
      end 
    end

    return invariants 
end

function findAllInvariants(domloops, uniqSet::Set{Symbol}, bbs)
  invariants = Dict{LivenessAnalysis.Loop, Set{Symbol} }()

  for l in domloops.loops
    invariants[l] = findInvariants(l.members, uniqSet)
  end

  return invariants
end

function processFuncCall(state :: memoizeState, func_expr, call_sig_arg_tuple)
  fetyp = typeof(func_expr)

  dprintln(3,"processFuncCall ", func_expr, " ", call_sig_arg_tuple, " ", fetyp)
  func = eval(func_expr)
  dprintln(3,"func = ", func, " type = ", typeof(func))

  ftyp = typeof(func)
  dprintln(4,"After name resolution: func = ", func, " type = ", ftyp)
  if ftyp == DataType
    return nothing
  end
  assert(ftyp == Function || ftyp == IntrinsicFunction || ftyp == LambdaStaticData)

  if ftyp == Function
    fs = (func, call_sig_arg_tuple)

    if haskey(state.mapNameFuncInfo, fs)
      return state.mapNameFuncInfo[fs]
    end

    dprintln(3,"Attempt to optimize ", func)
    ct = code_typed(func, call_sig_arg_tuple)
    dprintln(3,"ct = ", ct)
    dprintln(3,"typeof(ct) = ", typeof(ct))
    dprintln(3,"ct[1] = ", ct[1])
    dprintln(3,"typeof(ct[1]) = ", typeof(ct[1]))
    if typeof(ct[1]) == Expr
      ast = ct[1]
      dprintln(3,"ct head = ", ast.head) 
      assert(ast.head == :lambda)

      lives      = LivenessAnalysis.from_expr(ast)
      uniqSet    = AliasAnalysis.analyze_lambda(ast, lives)
      loop_info  = LivenessAnalysis.compute_dom_loops(lives)
#      invariants = findAllInvariants(loop_info, uniqSet, lives.basic_blocks)

      analyze_res = sparse_analyze(ast, lives, loop_info)
#      analyze_res = sparse_analyze(ast, lives, loop_info, invariants)
      if analyze_res != nothing
        state.mapNameFuncInfo[fs] = analyze_res
      end
      return analyze_res
    end
  elseif ftyp == LambdaStaticData
    fs = (func, call_sig_arg_tuple)

    if haskey(state.mapNameFuncInfo, fs)
      return state.mapNameFuncInfo[fs]
    end
    dprintln(3,"Adding ", func, " to functionsToProcess from LambdaStaticData")
    # FIX FIX FIX
    throw(string("Didn't add the case to handle LambdaStaticData in processFuncCall yet."))
  end
  return nothing
end

function opt_calls(x, state :: memoizeState, top_level_number, is_top_level, read)
  if typeof(x) == Expr
    if x.head == :call
      call_expr = x.args[1]
      call_sig = Expr(:tuple)
      call_sig.args = map(typeOfOpr, x.args[2:end])
      call_sig_arg_tuple = eval(call_sig)

      dprintln(2, "Start opt_calls = ", call_expr, " signature = ", call_sig_arg_tuple)

      for i = 2:length(x.args)
        new_arg = AstWalker.AstWalk(x.args[i], opt_calls, state)
        assert(isa(new_arg,Array))
        assert(length(new_arg) == 1)
        x.args[i] = new_arg[1]
      end

      process_res = processFuncCall(state, call_expr, call_sig_arg_tuple)

      if process_res != nothing
        x.args[1] = process_res
        println(2, "Replaced call_name = ", call_name, " type = ", typeof(call_name), " new = ", x.args[1])
      end

      return x
    end    
  end
  nothing
end

gSparseAccelerateState = memoizeState()

macro acc(ast)
  println(2, "Mtest ", ast, " ", typeof(ast), " gSparseAccelerateState = ", gSparseAccelerateState)
  res = AstWalker.AstWalk(ast, opt_calls, gSparseAccelerateState)
  assert(isa(res,Array))
  assert(length(res) == 1)
  return esc(res[1])
end

function foo(x)
  println("foo = ", x)
#  z1 = zeros(100)
#  sum = 0.0
#  for i = 1:100
#    sum = sum + z1[i] 
#  end
  x+1
end

function foo_new(x)
  println("foo_new = ", x)
  x+10
end

function bar(x)
  println("bar = ", x)

  z1 = zeros(100)
  sum = 0.0
  for i = 1:100
    sum = sum + z1[i] 
  end

  x+2
end

function bar_new(x)
  println("bar_new = ", x)
  x+20
end

y = 1
#@acc y = foo(bar(y))
#@acc y = foo(7)
println(y)
processFuncCall(gSparseAccelerateState, foo, (Int64,))

end   # end of module
