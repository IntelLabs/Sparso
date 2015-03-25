module SparseAccelerator

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

testing_mode = false
function enable_testing_mode()
  dprintln(1,"enabling testing mode")
  global testing_mode = true
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

function initSymbol2TypeDict(expr)
    assert(expr.head == :lambda) # (:lambda, {param, meta@{localvars, types, freevars}, body})

    local types = expr.args[2][2]
    symbolInfo = Dict{Symbol, Any}()
    for i = 1:length(types)
        symbolInfo[types[i][1]] = types[i][2]
        dprintln(4,"Add symbolInfo: ", types[i][1], "==>", symbolInfo[types[i][1]])
    end
    symbolInfo
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
    # use get() instead [] in case the key (like Symbol "*") does not exist
    return get(symbolInfo, ast, Nothing)
  else
    return asttyp
  end
end

# Our type inference annotates Expr or SymbolNode with the following types:
#   AbstractSparseMatrix, Matrix, Vector, Array, Bool, Number, Nothing

function checkDistributivityForCall(head, args, symbolInfo, distributive)
  println("\t**** checkDistributivityForCall", head, " ", args)
  if head == :(*)
    # "*" can have 2 or 3 operands. Example: [:*,:α,:A,:p]
    if length(args) == 3
        distributive = checkDistributivity(args[1], symbolInfo, distributive)
        return checkDistributivityForCall(head, args[2:end], symbolInfo, distributive)
    end 
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    if typeOfNode(args[1], symbolInfo) <: AbstractSparseMatrix
        if (typeOfNode(args[2], symbolInfo) <: AbstractSparseMatrix) || (typeOfNode(args[2], symbolInfo) <: Number)
            return distributive
        elseif (typeOfNode(args[2], symbolInfo) <: Vector)
            return distributive
        else
            return false
        end
    end
    if typeOfNode(args[1], symbolInfo) <: Number
        if (typeOfNode(args[2], symbolInfo) <: AbstractSparseMatrix)
            return distributive
        elseif (typeOfNode(args[2], symbolInfo) <: Vector)
            return distributive
        else
            return false
        end
    end
  end
  if head == :(+) || head == :(-) # Arithmetic operators
    if (length(args) == 1) 
        return checkDistributivity(args[1], symbolInfo, distributive)
    end
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    if typeOfNode(args[1], symbolInfo) <: AbstractSparseMatrix
        if (typeOfNode(args[2], symbolInfo) <: AbstractSparseMatrix) 
            return distributive
        else
            return false
        end
    end
    if typeOfNode(args[1], symbolInfo) <: Vector
        if (typeOfNode(args[2], symbolInfo) <: Vector)
            return distributive
        else 
            return false
        end
    end
  end
  if head == :(\)
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    if typeOfNode(args[1], symbolInfo) <: AbstractSparseMatrix
        if (typeOfNode(args[2], symbolInfo) <: Vector)
            return distributive
        else 
            return false
        end
    end
  end
  if head == :dot
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    if (typeOfNode(args[1], symbolInfo) <: Vector) && (typeOfNode(args[2], symbolInfo) <: Vector) 
        return distributive
    else
        return false
    end        
  end
  if head == :(==) || head == :(!=) || head == :(<) || head == :(<=) || 
     head == :(>) || head == :(>=) # Numeric comparison 
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    return distributive
  end
  if head == :(~) # Bitwise operator
    return checkDistributivity(args[1], symbolInfo, distributive)
  end
  if head == :(&) || head == :(|) || head == :($)
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    return distributive
  end
  if head == :(>>>) || head == :(>>) || head == :(<<) # Bitwise operators
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    return distributive
  end
  if head == :(+) || head == :(-) || head == :(*) || head == :(/) || head == :(\) ||
     head == :(^) || head == :(%) # Arithmetic operators
    # "+" and "-" have been tested before as unary operations. Here handle binary op case
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    return distributive
  end
  if head == :(.+) || head == :(.-) || head == :(.*) || head == :(./) || head == :(.\) ||
     head == :(.^) || head == :(.==) || head == :(.!=) || head == :(.<) ||
     head == :(.<=) || head == :(.>) || head == :(.>=) # Bitwise operators
    distributive = checkDistributivityForBinaryOP(args, symbolInfo, distributive)
    return distributive
  end
  if head == :(!) # Negation on bool
    return checkDistributivity(args[1], symbolInfo, distributive)
  end
  if head == :norm 
    return checkDistributivityForCall(:dot, [args[1], args[1]], symbolInfo, distributive)
  end
  if head == :sqrt
    return checkDistributivity(args[1], symbolInfo, distributive)   
  end
  if isa(head, TopNode)
    return checkDistributivityForCall(head.name, args, symbolInfo, distributive)
  end
  if head == :copy
    return checkDistributivity(args[1], symbolInfo, distributive)
  end
  if head == :diag
    return distributive
  end
  # This is to quickly pass pagerank. However, we need a more general machenism. we cannot
  # write all kinds of functions tediously.
  # TODO: a general mechanism to handle functions
   if head == :max || head == :vec || head == :sum || head == :colon || head == :start || head == :done || head == :next || head == :tupleref || head == :toc || head == :tic || head == :time || head == :clock_now
    return distributive
  end
 
  throw(string("checkDistributivityForCall: unknown AST (", head, ",", args, ")"))
end

function checkDistributivityForBinaryOP(ast, symbolInfo, distributive)
  assert(length(ast) == 2)
  local lhs = ast[1]
  local rhs = ast[2]
  distributive = checkDistributivity(rhs, symbolInfo, distributive)
  distributive = checkDistributivity(lhs, symbolInfo, distributive)
  distributive
end

function checkDistributivityForExprs(ast, symbolInfo, distributive)
  local len = length(ast)
  local typ = Nothing
  for i = 1:len
    distributive = checkDistributivity(ast[i], symbolInfo, distributive)
  end
  distributive
end

function checkDistributivityForIf(args, symbolInfo, distributive)
    # The structure of the if node is an array of length 2.
    # The first index is the conditional.
    # The second index is the label of the else section.
    assert(length(args) == 2)
    if_clause  = args[1]
    return checkDistributivity(if_clause, symbolInfo, distributive)
end

# :lambda expression
# ast = [ parameters, meta (local, types, etc), body ]
function checkDistributivityForLambda(ast, symbolInfo, distributive)
  assert(length(ast) == 3)
  local body  = ast[3]
  return checkDistributivity(body, symbolInfo, distributive)
end

# Recursively infer types of the nodes in the AST. Meanwhile, test for
# reordering distributivity      
# TODO: if only for reordering purpose, once we know distributive is false, we
# can stop inferring types, since we cannot do reordering optimization anyway.
function checkDistributivity(ast::Any, symbolInfo::Dict{Symbol, Any}, distributive::Bool)
  local asttyp = typeof(ast)
  dprintln(2,"checkDistributivity typeof=", asttyp, " typeOfNode=", typeOfNode(ast, symbolInfo), " ", ast)

  if isa(ast, Tuple)
    for i = 1:length(ast)
        distributive |= checkDistributivity(ast[i], symbolInfo, distributive)
    end
  elseif asttyp == Expr
    local head = ast.head
    local args = ast.args
    if head == :lambda
        distributive = checkDistributivityForLambda(args, symbolInfo, distributive)
        dprintln(2, "\tLambda ", distributive)
    elseif head == :body
        distributive = checkDistributivityForExprs(args, symbolInfo, distributive)
        dprintln(2, "\tBody ", distributive)
    elseif head == :(=) || head == :(+=) || head == :(-=) || head == :(*=) || 
           head == :(/=) || head == :(\=) || head == :(%=) || head == :(^=) ||
           head == :(&=) || head == :(|=) || head == :($=) || head == :(>>>=) ||
           head == :(>>=) || head == :(<<=)  # Updating operators
        assert(length(args) == 2)
        local lhs = args[1]
        local rhs = args[2]
        distributive = checkDistributivity(rhs, symbolInfo, distributive)
        dprintln(2, "\tAssignment ", distributive)
    elseif head == :return
        distributive = checkDistributivity(args, symbolInfo, distributive)
        dprintln(2, "\tReturn ", distributive)
    elseif head == :call || head == :call1
        distributive = checkDistributivityForCall(args[1], args[2:end], symbolInfo, distributive)
        dprintln(2, "\tCall ", distributive)
    elseif head == :gotoifnot
        distributive = checkDistributivityForIf(args, symbolInfo, distributive)
        dprintln(2, "\tGotoifnot ", distributive)
    elseif head == :line
    else
        throw(string("checkDistributivity: unknown Expr head :", head))
    end
    return distributive 
  elseif asttyp == SymbolNode
    dprintln(2, "\tSymbolNode ", distributive)
    return distributive
  elseif asttyp == Symbol
    dprintln(2, "\tSymbol ", distributive)
    return distributive
  elseif asttyp == LabelNode ||
         asttyp == GotoNode || asttyp == LineNumberNode ||
         asttyp == ASCIIString || asttyp == LambdaStaticData
  elseif asttyp <: Number
      dprintln(2, "\tNumber ", distributive)
      return distributive
  elseif asttyp.name == Array.name # Example: Array{Any,1}
      typ = asttyp.parameters[1]
      if (asttyp.parameters[2] == 1) #1-dimensonal array is vector
          typ = Vector
      elseif (asttyp.parameters[2] == 2) #2-dimensonal array is matrix
          typ = Matrix
      end
      dprintln(2, "\tArray ", distributive)
      return distributive
  else
      throw(string("checkDistributivity: unknown AST type :", asttyp, " ", ast))
  end
  dprintln(2, "\tOther ", distributive)
  return distributive
end

function SparseOptimize(ast, call_sig_arg_tuple, call_sig_args)
  assert(typeof(ast)== Expr)
  assert(ast.head == :lambda)

  dprintln(3,"SparseOptimize args = ", call_sig_arg_tuple, "\n", ast, "\n")

  # Populate the lambda's meta info for symbols with what we know from the typed AST
  # This is necessary, since the AST may contain symbols, which does not contain
  # type info itself. Example: x= ... That x is a Symbol, not a SymbolNode, and
  # thus does not have type info stored. We have to look up from the lambda. To be
  # faster, we build this dictionary, and look up from it instead.
  symbolInfo = initSymbol2TypeDict(ast)

  distributive = checkDistributivity(ast, symbolInfo, true)
  dprintln(3,"After our type inference, distributive = ", distributive)
      
  if !distributive
    return ast
  end

  body = ast.args[3]

  lives = LivenessAnalysis.initial_expr(ast)
  dprintln(3,"function to analyze type = ", typeof(body.args), "\n", body)
  dprintln(3,"testing_mode = ", testing_mode)
  if testing_mode
    LivenessAnalysis.insertBetween(lives, 2, 4)
    #LivenessAnalysis.insertBefore(lives, 2)
    body_reconstructed = LivenessAnalysis.createFunctionBody(lives)
    dprintln(3,"reconstructed_body type = ", typeof(body_reconstructed.args), "\n", body_reconstructed)
  end
  loop_info  = LivenessAnalysis.compute_dom_loops(lives)

  analyze_res = sparse_analyze(ast, lives, loop_info, symbolInfo)
  dprintln(3,"result after sparse_analyze\n", analyze_res, " type = ", typeof(analyze_res))
  assert(typeof(analyze_res) == Expr)
  assert(analyze_res.head == :lambda)
  return analyze_res
end

LivenessAnalysis.set_debug_level(4)

end   # end of module
