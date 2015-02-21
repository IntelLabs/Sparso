module SparseAccelerator

export @acc

include("ast_walk.jl")
include("liveness.jl")
include("alias-analysis.jl")

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
#  dprintln(3,"typeOfOpr ", x, " type = ", typeof(x))
  if isa(x, Expr) x.typ
  elseif isa(x, SymbolNode) x.typ
  else typeof(x) 
  end
end   

type memoizeState
  mapNameFuncInfo :: Dict{Any, Any}   # tracks the mapping from unoptimized function name to optimized function name
  trampolineSet   :: Set{Any}         # tracks whether we've previously created a trampoline for a given function name and signature

  function memoizeState()
    new (Dict{Any,Any}(), Set{Any}())
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
    invariants[l] = findInvariants(l.members, uniqSet, bbs)
  end

  return invariants
end

function buildSymbol2TypeDict(expr)
    assert(expr.head == :lambda) # (:lambda, {param, meta@{localvars, types, freevars}, body})
    local meta  = expr.args[2]
    dprintln(3, "meta is:", meta)
    symbolInfo = Dict{Symbol, Any}()
    for info in meta[2]
        symbolInfo[info[1]] = info
        dprintln(3,"Add symbolInfo: ", info[1], "==>", symbolInfo[info[1]])
    end
    dprintln(3, "SI is ", symbolInfo)
    symbolInfo
end

function updateLambdaMeta(expr, symbolInfo)
    assert(expr.head == :lambda) # (:lambda, {param, meta@{localvars, types, freevars}, body})
    
    # First, populate the type info of symbols into the metadata.
    local meta  = expr.args[2]
    dprintln(3, "meta is:", meta)
    
    for i in 1:length(meta[2])
        info = meta[2][i]
        #TODO: graceful exit of analysis in case the key is not in the dictionary
        meta[2][i] = symbolInfo[info[1]]
    end
    expr.args[2] = meta
    dprintln(3, "Lambda is updated:", expr)
end

# In our type inference, we put inferred types into the typ field for expressions.
# A symbol's type is still looked up from a dictionary.
# TODO for Todd: when generating a new AST from CFG after our optimizations, ignore
# those fields set up by our type inference.
function typeOfNode(ast, symbolInfo) 
  local asttyp = typeof(ast)
  if asttyp == Expr
    return ast.typ
  elseif asttyp == SymbolNode
    return ast.typ
  elseif asttyp == Symbol
    return symbolInfo[ast][2]
  else 
    return nothing
  end
end
    
function inferTypesForCall(head, args, symbolInfo, distributive)
  if head == :(*)
    if typeOfNode(args[1], symbolInfo) <: AbstractSparseMatrix
        if (typeOfNode(args[2], symbolInfo) <: AbstractSparseMatrix) ||
           (typeOfNode(args[2], symbolInfo) <: Number)
            return (AbstractSparseMatrix, distributive)
        elseif (typeOfNode(args[2], symbolInfo) <: Vector)
            return (Vector, distributive)
        else
            return (Any, false)
        end
    end
    if typeOfNode(args[1], symbolInfo) <: Number
        if (typeOfNode(args[2], symbolInfo) <: AbstractSparseMatrix)
            return (AbstractSparseMatrix, distributive)
        elseif (typeOfNode(args[2], symbolInfo) <: Vector)
            return (Vector, distributive)
        else
            return (Any, false)
        end
    end
  end
  if head == :(+) || head == :(-)
    if typeOfNode(args[1], symbolInfo) <: AbstractSparseMatrix
        if (typeOfNode(args[2], symbolInfo) <: AbstractSparseMatrix) 
            return (AbstractSparseMatrix, distributive)
        else
            return (Any, false)
        end
    end
    if typeOfNode(args[1], symbolInfo) <: Vector
        if (typeOfNode(args[2], symbolInfo) <: Vector)
            return (Vector, distributive)
        else 
            return (Any, false)
        end
    end
  end
  if head == :(\)
    if typeOfNode(args[1], symbolInfo) <: AbstractSparseMatrix
        if (typeOfNode(args[2], symbolInfo) <: Vector)
            return (Vector, distributive)
        else 
            return (Any, false)
        end
    end
  end
  if head == :dot
    if (typeOfNode(args[1], symbolInfo) <: Vector) && (typeOfNode(args[2], symbolInfo) <: Vector) 
        return (Number, distributive)
    else
        return (Any, false)
    end        
  end
  if head == :(==) || head == :(!=) || head == :(<) || head == :(<=) || 
     head == :(>) || head == :(>=) # Numeric comparison 
    return inferTypesForBinaryOP(args, symbolInfo, distributive)
  end
  if head == :(~) # Bitwise operator
    return inferTypes(args[1], symbolInfo, distributive)
  end
  if head == :(&) || head == :(|) || head == :($) ||
     head == :(>>>) || head == :(>>) || head == :(<<) # Bitwise operators
    return inferTypesForBinaryOP(args, symbolInfo, distributive)
  end
  if head == :(+) || head == :(-) # Arithmetic operators
    if (length(args) == 1) 
        return inferTypes(args[1], symbolInfo, distributive)
    else
        return inferTypesForBinaryOP(args, symbolInfo, distributive)
    end
  end
  if head == :(*) || head == :(/) || head == :(\) ||  # Arithmetic operators
     head == :(^) || head == :(%) 
    return inferTypesForBinaryOP(args, symbolInfo, distributive)
  end
  if head == :(!) # Negation on bool
    return inferTypes(args[1], symbolInfo, distributive)
  end
  if head == :norm 
    return inferTypes(args, symbolInfo, distributive)
  end
  if isa(head, TopNode)
    return inferTypesForCall(head.name, args, symbolInfo, distributive)
  end
  throw(string("inferTypesForCall: unknown AST (", head, ",", args, ")"))
end

function inferTypesForBinaryOP(ast, symbolInfo, distributive)
  assert(length(ast) == 2)
  local lhs = ast[1]
  local rhs = ast[2]
  (typ, distributive) = inferTypes(rhs, symbolInfo, distributive)
  (typ, distributive) = inferTypes(lhs, symbolInfo, distributive)
  (typ, distributive)
end

function inferTypesForExprs(ast, symbolInfo, distributive)
  local len = length(ast)
  local typ
  for i = 1:len
    (typ, distributive) = inferTypes(ast[i], symbolInfo, distributive)
  end
  (typ, distributive)
end

function inferTypesForIf(args, symbolInfo, distributive)
    # The structure of the if node is an array of length 2.
    # The first index is the conditional.
    # The second index is the label of the else section.
    assert(length(args) == 2)
    if_clause  = args[1]

    return inferTypes(if_clause, symbolInfo, distributive)
end

# :lambda expression
# ast = [ parameters, meta (local, types, etc), body ]
function inferTypesForLambda(ast, symbolInfo, distributive)
  assert(length(ast) == 3)
  local body  = ast[3]
  return inferTypes(body, symbolInfo, distributive)
end

# Recursively infer types of the nodes in the AST. Meanwhile, test for
# reordering distributivity      
# TODO: if only for reordering purpose, once we know distributive is false, we
# can stop inferring types, since we cannot do reordering optimization anyway.
function inferTypes(ast::Any, symbolInfo::Dict{Symbol, Any}, distributive::Bool)
  local asttyp = typeof(ast)
  dprintln(2,"inferTypes ", asttyp)

  if isa(ast, Tuple)
    for i = 1:length(ast)
        distributive |= inferTypes(ast[i], symbolInfo, distributive)
    end
  elseif asttyp == Expr
    local head = ast.head
    local args = ast.args
    dprintln(2, "Expr ", head, " ", args)
    if head == :lambda
        return inferTypesForLambda(args, symbolInfo, distributive)
    elseif head == :body
        return inferTypesForExprs(args, symbolInfo, distributive)
    elseif head == :(=) || head == :(+=) || head == :(-=) || head == :(*=) || 
           head == :(/=) || head == :(\=) || head == :(%=) || head == :(^=) ||
           head == :(&=) || head == :(|=) || head == :($=) || head == :(>>>=) ||
           head == :(>>=) || head == :(<<=)  # Updating operators
        return inferTypesForBinaryOP(args, symbolInfo, distributive)
    elseif head == :return
        return inferTypes(args, symbolInfo, distributive)
    elseif head == :call || head == :call1
        return inferTypesForCall(args[1], args[2:end], symbolInfo, distributive)
    elseif head == :gotoifnot
        return inferTypesForIf(args, symbolInfo, distributive)
    elseif head == :line
    else
        throw(string("inferTypes: unknown Expr head :", head))
    end
  elseif asttyp == SymbolNode
    ast.typ = symbolInfo[ast.Name][2]
    return (ast.typ, distributive)
  elseif asttyp == Symbol || asttyp == SymbolNode || asttyp == LabelNode ||
         asttyp == GotoNode || asttyp == LineNumberNode ||
         asttyp == ASCIIString || asttyp == LambdaStaticData
  elseif asttyp <: Number
      return (Number, distributive)
  elseif asttyp.name == Array.name # Example: Array{Any,1}
      return (Vector, distributive)
  else
      throw(string("inferTypes: unknown AST type :", asttyp, " ", ast))
  end
  return (asttyp, distributive)
end

function processFuncCall(func_expr, call_sig_arg_tuple)
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

    dprintln(3,"Attempt to optimize ", func)
    ct = code_typed(func, call_sig_arg_tuple)
    dprintln(3,"ct = ", ct)
    dprintln(3,"typeof(ct) = ", typeof(ct))
    dprintln(3,"ct[1] = ", ct[1])
    dprintln(3,"typeof(ct[1]) = ", typeof(ct[1]))
    if typeof(ct[1]) == Expr
      # We get the type info from the typed AST, but do analysis on the lowered AST.
      # Lowered AST is close to math level, while typed AST may contain detailed
      # internal implementations of the math expressions due to optimizations done
      # during type inference (such as inlining, tuple elimination, etc.)
      # Our analyses target math expressions, not depend on how they are implemented.
      symbolInfo = buildSymbol2TypeDict(ct[1])

      cl = code_lowered(func, call_sig_arg_tuple)
      dprintln(3,"cl = ", cl)
      dprintln(3,"typeof(cl) = ", typeof(cl))
      dprintln(3,"cl[1] = ", cl[1])
      dprintln(3,"typeof(cl[1]) = ", typeof(cl[1]))
           
      ast = cl[1]
      dprintln(3,"cl head = ", ast.head) 
      assert(ast.head == :lambda)

      # Populate the lambda's meta info for symbols with what we know from the typed AST
      updateLambdaMeta(ast, symbolInfo)

      # Unfortunately, the nodes of typed AST and lowered AST may not be 1:1 correspondence.
      # Although we know the types of the symbols from typed AST, we do not know the types 
      # of the (intermediate) nodes of lowered AST. So do a simple type inference on the 
      # lowered AST
      (typ, distributive) = inferTypes(ast.args[3], symbolInfo, true)
      if !distributive
        return nothing
      end

      lives      = LivenessAnalysis.from_expr(ast)
#      uniqSet    = AliasAnalysis.analyze_lambda(ast, lives)
      loop_info  = LivenessAnalysis.compute_dom_loops(lives)
#      invariants = findAllInvariants(loop_info, uniqSet, lives.basic_blocks)

      analyze_res = sparse_analyze(ast, lives, loop_info)
#      analyze_res = sparse_analyze(ast, lives, loop_info, invariants)
      return analyze_res
    end
  elseif ftyp == LambdaStaticData
    fs = (func, call_sig_arg_tuple)

    dprintln(3,"Adding ", func, " to functionsToProcess from LambdaStaticData")
    # FIX FIX FIX
    throw(string("Didn't add the case to handle LambdaStaticData in processFuncCall yet."))
  end
  return nothing
end

gSparseAccelerateState = memoizeState()

function opt_calls_insert_trampoline(x, state :: memoizeState, top_level_number, is_top_level, read)
  if typeof(x) == Expr
    if x.head == :call
      call_expr = x.args[1]
      call_sig_args = x.args[2:end]
      dprintln(2, "Start opt_calls = ", call_expr, " signature = ", call_sig_args, " typeof(call_expr) = ", typeof(call_expr))

      new_func_name = string("opt_calls_trampoline_", string(call_expr))
      new_func_sym  = symbol(new_func_name)

      for i = 2:length(x.args)
        new_arg = AstWalker.AstWalk(x.args[i], opt_calls_insert_trampoline, state)
        assert(isa(new_arg,Array))
        assert(length(new_arg) == 1)
        x.args[i] = new_arg[1]
      end

      tmtup = (call_expr, call_sig_args)
      if !in(tmtup, state.trampolineSet)
        dprintln(3,"Creating new trampoline for ", call_expr)
        push!(state.trampolineSet, tmtup)
        println(new_func_sym)
        for i = 1:length(call_sig_args)
          println("    ", call_sig_args[i])
        end
 
       @eval function ($new_func_sym)(orig_func, $(call_sig_args...))
              call_sig = Expr(:tuple)

              call_sig.args = map(typeOfOpr, Any[ $(call_sig_args...) ]) 
              call_sig_arg_tuple = eval(call_sig)
              println(call_sig_arg_tuple)

              fs = ($new_func_sym, call_sig_arg_tuple)

              if haskey(gSparseAccelerateState.mapNameFuncInfo, fs)
                func_to_call = gSparseAccelerateState.mapNameFuncInfo[fs]
              else
                process_res = processFuncCall(orig_func, call_sig_arg_tuple)

                if process_res != nothing
                  dprintln(3,"processFuncCall DID optimize ", orig_func)
                  func_to_call = process_res
                else
                  dprintln(3,"processFuncCall didn't optimize ", orig_func)
                  func_to_call = orig_func
                end
                gSparseAccelerateState.mapNameFuncInfo[fs] = func_to_call
              end

              println("running ", $new_func_name, " fs = ", fs)
              func_to_call($(call_sig_args...))
            end
      end

      resolved_name = @eval SparseAccelerator.$new_func_sym
      x.args = [ resolved_name, call_expr, x.args[2:end] ]

      dprintln(2, "Replaced call_expr = ", call_expr, " type = ", typeof(call_expr), " new = ", x.args[1])

      return x
    end    
  end
  nothing
end


function convert_expr(ast)
  dprintln(2, "Mtest ", ast, " ", typeof(ast), " gSparseAccelerateState = ", gSparseAccelerateState)
  res = AstWalker.AstWalk(ast, opt_calls_insert_trampoline, gSparseAccelerateState)
  assert(isa(res,Array))
  assert(length(res) == 1)
  dprintln(2,res[1])
  return esc(res[1])
end


macro acc(ast)
  convert_expr(ast)
end


function foo(x)
  println("foo = ", x, " type = ", typeof(x))
  z1 = zeros(100)
  sum = 0.0
  for i = 1:100
    sum = sum + z1[i] 
  end
#  bt = backtrace()
#  Base.show_backtrace(STDOUT, bt)
#  println("")
  x+1
end

function foo2(x,y)
  println("foo2")
  x+y
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

function testit()
#y = 1
z = 7
#@acc y = foo(bar(y))
#println(macroexpand(quote @acc y = foo(z) end))
#@acc y = foo(z)
#@acc y = foo(y)
println(y)
#processFuncCall(foo, (Int64,))
#convert_expr(quote y = foo(z) end)
end

#testit()

end   # end of module
