module SparseAccelerator

using CompilerTools

include("alias-analysis.jl")

# This controls the debug print level.  0 prints nothing.  At the moment, 2 prints everything.
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
  invariants = Dict{CompilerTools.LivenessAnalysis.Loop, Set{Symbol} }()

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
#  println("\t**** checkDistributivityForCall ", head, " ", args)
  if head == :(*) || head == :SpMV
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
  if head == :Array
    return distributive
  end
  if isa(head, Expr)
    dprintln(1, "head.head = ", head.head)
    if head.head == :call
      dprintln(1, "call")
      dprintln(1, head.args[1], " type = ", typeof(head.args[1]))
      if isa(head.args[1], TopNode)
        dprintln(1, "TopNode ", head.args[1:end])
        if head.args[2] == :SparseAccelerator
          return distributive
        end 
      end
    end
  end
#  if head == Top
#    println("getfield ", args[1], " ", args[2], " ", typeof(args[1]), " " , typeof(args[2]))
#    #return distributive
#  end
  # This is to quickly pass pagerank. However, we need a more general machenism. we cannot
  # write all kinds of functions tediously.
  # TODO: a general mechanism to handle functions
   if head == :max || head == :vec || head == :sum || head == :colon || head == :start || head == :done || head == :next || head == :tupleref || head == :toc || head == :tic || head == :time || head == :clock_now || head == :println ||
      head == :spones ||  # spones,Any[:(A::Union((Int64,Int64,Int64,Any,Any,Any),Base.SparseMatrix.SparseMatrixCSC{Float64,Int32}))]
      head == :size || # size,Any[:(A::Base.SparseMatrix.SparseMatrixCSC{Float64,Int32}),1]
      head == :repmat || # repmat,Any[:((top(vect))(1 / m::Int64::Float64)::Array{Float64,1}),:(m::Int64)]
      head == :scale || # scale,Any[:(A::Base.SparseMatrix.SparseMatrixCSC{Float64,Int32}),:(1 ./ d::Array{Float64,1}::Array{Float64,1})]
      head == :getfield
    return distributive
  end

  println("head type = ", typeof(head)) 
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
  elseif asttyp == NewvarNode
    dprintln(2, "\tNewvarNode ", distributive)
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

  body = ast.args[3]

  lives = CompilerTools.LivenessAnalysis.from_expr(ast)
  dprintln(3,"function to analyze type = ", typeof(body.args), "\n", body)
  dprintln(3,"testing_mode = ", testing_mode)
  if testing_mode
    CompilerTools.LivenessAnalysis.insertBetween(lives, 2, 4)
    #CompilerTools.LivenessAnalysis.insertBefore(lives, 2)
    body_reconstructed = Expr(:body)
    body_reconstructed.typ = body.typ
    body_reconstructed.args = CompilerTools.LivenessAnalysis.createFunctionBody(lives)
    dprintln(3,"reconstructed_body type = ", typeof(body_reconstructed.args), "\n", body_reconstructed)
  end
  loop_info  = CompilerTools.Loops.compute_dom_loops(lives)

  analyze_res = sparse_analyze(ast, lives, loop_info, symbolInfo)
  dprintln(3,"result after sparse_analyze\n", analyze_res, " type = ", typeof(analyze_res))
  assert(typeof(analyze_res) == Expr && analyze_res.head == :lambda)
  dprintln(3,"typeof(args[3]) = ", typeof(analyze_res.args[3]))
  if typeof(analyze_res.args[3]) == Expr
    dprintln(3,"args[3].head = ", analyze_res.args[3].head)
  end
  assert(typeof(analyze_res.args[3]) == Expr && analyze_res.args[3].head == :body)
  return analyze_res
end

CompilerTools.LivenessAnalysis.set_debug_level(0)

# Choose between Julia and PCL library. 
const JULIA_LIB = 0
const PCL_LIB = 1
DEFAULT_LIBRARY = JULIA_LIB
function use_lib(value) 
  global DEFAULT_LIBRARY = value
end

function CreateCSR(A::SparseMatrixCSC)
  ccall((:CSR_Create, LIB_PATH), Ptr{Void},
       (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint),
       A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval), 1)
end

function DestroyCSR(A::Ptr{Void})
  ccall((:CSR_Destroy, LIB_PATH), Void,
        (Ptr{Void},),
        A)
end

# w = alpha*A*x + beta*y + gamma
function SpMV!(w::AbstractVector, alpha::Number, A::SparseMatrixCSC, x::AbstractVector, beta::Number, y::AbstractVector, gamma::Number)
    assert(length(w) == length(y))

    if DEFAULT_LIBRARY == PCL_LIB
        A2 = CreateCSR(A)
        ccall((:CSR_MultiplyWithVector, LIB_PATH), Void,
              (Ptr{Cdouble}, Cdouble, Ptr{Void}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
              pointer(w), alpha, A2, pointer(x), beta, pointer(y), gamma)
        DestroyCSR(A2)
    else
        # use Julia implementation
        w = alpha * A * x + beta * x + gamma
    end
end

# y = A*x
SpMV!(y::AbstractVector, A::SparseMatrixCSC, x::AbstractVector) = SpMV!(y, one(eltype(x)), A, x, zero(eltype(y)), y, zero(eltype(x)))

function PageRank!(w::AbstractVector, alpha::Number, A::SparseMatrixCSC, x::AbstractVector, beta::Number, y::AbstractVector, gamma::Number, z::AbstractVector)
    assert(length(w) == length(y))
    assert(length(w) == length(z))

    if DEFAULT_LIBRARY == PCL_LIB
        A2 = CreateCSR(A)
        ccall((:PageRank, LIB_PATH), Void,
              (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
              length(w), pointer(w), alpha, pointer(A.colptr), pointer(A.rowval), pointer(x), beta, pointer(y), gamma, pointer(z))
        DestroyCSR(A2)
    else
        # use Julia implementation
        w = (alpha * A * x + beta * x + gamma).*z
    end
end

# alpha*A*x + beta*y + gamma
function SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector, beta::Number, y::AbstractVector, gamma::Number)
  w = Array(Cdouble, length(x))
  SpMV!(w, alpha, A, x, beta, y, gamma)
  w
end

# A*x
SpMV(A::SparseMatrixCSC, x::AbstractVector) = SpMV(one(eltype(x)), A, x, zero(eltype(x)), x, zero(eltype(x)))

function WAXPBY!(w::Vector, alpha::Number, x::Vector, beta::Number, y::Vector)
  assert(length(x) == length(y))
  assert(length(x) == length(w))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:waxpby, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
          length(x), pointer(w), alpha, pointer(x), beta, pointer(y))
  else
    w = alpha*x + beta*y
  end
end

function WAXPBY(alpha::Number, x::Vector, beta::Number, y::Vector)
  w = Array(Cdouble, length(x))
  WAXPBY!(w, alpha, x, beta, y)
  w
end

function Dot(x::Vector, y::Vector)
  assert(length(x) == length(y))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:dot, LIB_PATH), Cdouble,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(x), pointer(y))
  else
    dot(x, y)
  end
end

function PointwiseDivide!(w::Vector, x::Vector, y::Vector)
  assert(length(x) == length(y))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:pointwiseDivide, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(w), pointer(x), pointer(y))
  else
    w = x./y
  end
end

function PointwiseDivide(x::Vector, y::Vector)
  w = Array(Cdouble, length(x))
  PointwiseDivide!(w, x, y)
  w
end

function PointwiseMultiply!(w::Vector, x::Vector, y::Vector)
  assert(length(x) == length(y))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:pointwiseMultiply, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          length(x), pointer(w), pointer(x), pointer(y))
  else
    w = x.*y
  end
end

function PointwiseMultiply(x::Vector, y::Vector)
  w = Array(Cdouble, length(x))
  PointwiseMultiply!(w, x, y)
  w
end

function WAXPB!(w::Vector, alpha::Number, x::Vector, beta::Number)
  assert(length(w) == length(x))

  if DEFAULT_LIBRARY == PCL_LIB
    ccall((:waxpb, LIB_PATH), Void,
          (Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
          length(x), pointer(w), alpha, pointer(x), beta)
  else
    w = alpha*x + beta
  end
end

function WAXPB(alpha::Number, x::Vector, beta::Number)
  w = Array(Cdouble, length(x))
  WAXPB!(w, alpha, x, beta)
  w
end

end   # end of module

#function Base.A_mul_B!(alpha::Number, A::SparseMatrixCSC, x::AbstractVector, beta::Number, y::AbstractVector)
#  SparseAccelerator.SpMV(alpha, A, x, beta, y)
#end

