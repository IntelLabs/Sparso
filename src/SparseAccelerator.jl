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

function typeOfOpr(x)
#  dprintln(3,"typeOfOpr ", x, " type = ", typeof(x))
  if isa(x, Expr) x.typ
  elseif isa(x, SymbolNode) x.typ
  else typeof(x) 
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

function initSymbol2TypeDict(expr, call_sig_arg_tuple)
    assert(expr.head == :lambda) # (:lambda, {param, meta@{localvars, types, freevars}, body})
    
    local param = expr.args[1]
    symbolInfo = Dict{Symbol, Any}()
    for i = 1:length(param)
        symbolInfo[param[i]] = call_sig_arg_tuple[i]
        dprintln(4,"Add symbolInfo: ", param[i], "==>", symbolInfo[param[i]])
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
    return symbolInfo[ast]
  else
    return asttyp
  end
end

function assignType(lhs, rhs, symbolInfo)
    # Usually, the left hand side should be a symbol(node) or a tuple of symbol(node)s
    if isa(lhs, Tuple)
        for i = 1:length(lhs)
                assignType(lhs[i], typeOfNode(rhs[i], symbolInfo), symbolInfo)
        end
    elseif typeof(lhs) == SymbolNode
        ast.typ = typeOfNode(rhs, symbolInfo)
    elseif typeof(lhs) == Symbol
        symbolInfo[lhs] = typeOfNode(rhs, symbolInfo)
    else
        throw("Unhandled LHS in assignment ", lhs, " type is ", typeof(lhs))
    end
    nothing
end

# Our type inference annotates Expr or SymbolNode with the following types:
#   AbstractSparseMatrix, Matrix, Vector, Array, Bool, Number, Nothing

function inferTypesForCall(head, args, symbolInfo, distributive)
  if head == :(*)
    # "*" can have 2 or 3 operands. Example: [:*,:α,:A,:p]
    if length(args) == 3
        (typ, distributive) = inferTypes(args[1], symbolInfo, distributive)
        assert(typ <: Number)
        return inferTypesForCall(head, args[2:end], symbolInfo, distributive)
    end 
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    if ltyp <: AbstractSparseMatrix
        if (rtyp <: AbstractSparseMatrix) || (rtyp <: Number)
            return (AbstractSparseMatrix, distributive)
        elseif (rtyp <: Vector)
            return (Vector, distributive)
        else
            return (Nothing, false)
        end
    end
    if ltyp <: Number
        if (rtyp <: AbstractSparseMatrix)
            return (AbstractSparseMatrix, distributive)
        elseif (rtyp <: Vector)
            return (Vector, distributive)
        else
            return (Nothing, false)
        end
    end
  end
  if head == :(+) || head == :(-) # Arithmetic operators
    if (length(args) == 1) 
        return inferTypes(args[1], symbolInfo, distributive)
    end
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    if ltyp <: AbstractSparseMatrix
        if (rtyp <: AbstractSparseMatrix) 
            return (AbstractSparseMatrix, distributive)
        else
            return (Nothing, false)
        end
    end
    if ltyp <: Vector
        if (rtyp <: Vector)
            return (Vector, distributive)
        else 
            return (Nothing, false)
        end
    end
  end
  if head == :(\)
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    if ltyp <: AbstractSparseMatrix
        if (rtyp <: Vector)
            return (Vector, distributive)
        else 
            return (Nothing, false)
        end
    end
  end
  if head == :dot
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    if (ltyp <: Vector) && (rtyp <: Vector) 
        return (Number, distributive)
    else
        return (Nothing, false)
    end        
  end
  if head == :(==) || head == :(!=) || head == :(<) || head == :(<=) || 
     head == :(>) || head == :(>=) # Numeric comparison 
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    return (Bool, distributive)
  end
  if head == :(~) # Bitwise operator
    return inferTypes(args[1], symbolInfo, distributive)
  end
  if head == :(&) || head == :(|) || head == :($)
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    return (ltyp, distributive)
  end
  if head == :(>>>) || head == :(>>) || head == :(<<) # Bitwise operators
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    return (ltyp, distributive)
  end
  if head == :(+) || head == :(-) || head == :(*) || head == :(/) || head == :(\) ||
     head == :(^) || head == :(%) # Arithmetic operators
    # "+" and "-" have been tested before as unary operations. Here handle binary op case
    (ltyp, rtyp, distributive) = inferTypesForBinaryOP(args, symbolInfo, distributive)
    return (ltyp, distributive)
  end
  if head == :(!) # Negation on bool
    return inferTypes(args[1], symbolInfo, distributive)
  end
  if head == :norm 
    return inferTypesForCall(:dot, [args[1], args[1]], symbolInfo, distributive)
  end
  if head == :sqrt
    return inferTypes(args[1], symbolInfo, distributive)   
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
  (rtyp, distributive) = inferTypes(rhs, symbolInfo, distributive)
  (ltyp, distributive) = inferTypes(lhs, symbolInfo, distributive)
  (ltyp, rtyp, distributive)
end

function inferTypesForExprs(ast, symbolInfo, distributive)
  local len = length(ast)
  local typ = Nothing
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
  dprintln(2,"inferTypes ", asttyp, " ", ast)

  if isa(ast, Tuple)
    for i = 1:length(ast)
        distributive |= inferTypes(ast[i], symbolInfo, distributive)
    end
  elseif asttyp == Expr
    local head = ast.head
    local args = ast.args
    local typ = Nothing
    if head == :lambda
        (typ, distributive) = inferTypesForLambda(args, symbolInfo, distributive)
        dprintln(2, "\tLambda (", typ, ", ", distributive, ")")
    elseif head == :body
        (typ, distributive) = inferTypesForExprs(args, symbolInfo, distributive)
        dprintln(2, "\tBody (", typ, ", ", distributive, ")")
    elseif head == :(=) || head == :(+=) || head == :(-=) || head == :(*=) || 
           head == :(/=) || head == :(\=) || head == :(%=) || head == :(^=) ||
           head == :(&=) || head == :(|=) || head == :($=) || head == :(>>>=) ||
           head == :(>>=) || head == :(<<=)  # Updating operators
        assert(length(args) == 2)
        local lhs = args[1]
        local rhs = args[2]
        (rtyp, distributive) = inferTypes(rhs, symbolInfo, distributive)
        assignType(lhs, rhs, symbolInfo)
        typ = rtyp
        dprintln(2, "\tAssignment (", typ, ", ", distributive, ")")
    elseif head == :return
        (typ, distributive) = inferTypes(args, symbolInfo, distributive)
        dprintln(2, "\tReturn (", typ, ", ", distributive, ")")
    elseif head == :call || head == :call1
        (typ, distributive) = inferTypesForCall(args[1], args[2:end], symbolInfo, distributive)
        dprintln(2, "\tCall (", typ, ", ", distributive, ")")
    elseif head == :gotoifnot
        (typ, distributive) = inferTypesForIf(args, symbolInfo, distributive)
        dprintln(2, "\tGotoifnot (", typ, ", ", distributive, ")")
    elseif head == :line
    else
        throw(string("inferTypes: unknown Expr head :", head))
    end
    ast.typ = typ
    return (typ, distributive) 
  elseif asttyp == SymbolNode
    ast.typ = symbolInfo[ast.Name]
    dprintln(2, "\tSymbolNode (", ast.typ, ", ", distributive, ")")
    return (ast.typ, distributive)
  elseif asttyp == Symbol
    typ = symbolInfo[ast]
    dprintln(2, "\tSymbol (", typ, ", ", distributive, ")")
    return (typ, distributive)
  elseif asttyp == LabelNode ||
         asttyp == GotoNode || asttyp == LineNumberNode ||
         asttyp == ASCIIString || asttyp == LambdaStaticData
  elseif asttyp <: Number
      dprintln(2, "\tNumber (Number, ", distributive, ")")
      return (Number, distributive)
  elseif asttyp.name == Array.name # Example: Array{Any,1}
      typ = asttyp.parameters[1]
      if (asttyp.parameters[2] == 1) #1-dimensonal array is vector
          typ = Vector
      elseif (asttyp.parameters[2] == 2) #2-dimensonal array is matrix
          typ = Matrix
      end
      dprintln(2, "\tArray (", typ, ", ", distributive, ")")
      return (typ, distributive)
  else
      throw(string("inferTypes: unknown AST type :", asttyp, " ", ast))
  end
  dprintln(2, "\tOther (", asttyp, ", ", distributive, ")")
  return (asttyp, distributive)
end

function SparseOptimize(ast, call_sig_arg_tuple)
  dprintln(3,"SparseOptimize args = ", call_sig_arg_tuple, "\n", ast, "\n")

  # Populate the lambda's meta info for symbols with what we know from the typed AST
  symbolInfo = initSymbol2TypeDict(ast, call_sig_arg_tuple)

  # Unfortunately, the nodes of typed AST and lowered AST may not be 1:1 correspondence.
  # Although we know the types of the symbols from typed AST, we do not know the types 
  # of the (intermediate) nodes of lowered AST. So do a simple type inference on the 
  # lowered AST
  (typ, distributive) = inferTypes(ast, symbolInfo, true)
  dprintln(3,"After our type inference, typ = ", typ, " distributive = ", distributive)
      
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
  return analyze_res
end

LivenessAnalysis.set_debug_level(4)

end   # end of module
