module SparseAccelerator

using CompilerTools

include("alias-analysis.jl")
include("function-description.jl")

# This controls the debug print level.  0 prints nothing.  At the moment, 2 prints everything.
DEBUG_LVL=0

# If true, keep analyzing, transformation, and execution as far as possible, 
# even when there is an error. This is to expose as many issues as possible, so
# that we can fix them collectively.
KEEP_GOING = false

# if true, for a loop region containing no other more expensive sparse matrix operations
# than SpMV, turn on reordering only when spmv is likely to be sped up by reordering 
USE_SPMV_REORDERING_POTENTIAL_MODEL = false

function set_debug_level(x)
    global DEBUG_LVL = x
end

function set_keep_going(x :: Bool)
    global KEEP_GOING = x
end

function set_use_spmv_reordering_potential_model(x :: Bool)
    global USE_SPMV_REORDERING_POTENTIAL_MODEL = x
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

# Create a path to libcsr.
const libcsr = joinpath(dirname(@__FILE__), "..", "lib", "libcsr.so")

include("sparse-analyze.jl")
include("exceptions.jl")

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

function is_assignment(expr :: Expr)
    head = expr.head
    return head == :(=) || head == :(+=) || head == :(-=) || head == :(*=) || 
           head == :(/=) || head == :(\=) || head == :(%=) || head == :(^=) ||
           head == :(&=) || head == :(|=) || head == :($=) || head == :(>>>=) ||
           head == :(>>=) || head == :(<<=)
end

function memoize_GenSym_types(ast, symbolInfo, top_level_number, is_top_level, read)
    if typeof(ast) == Expr && is_assignment(ast) && typeof(ast.args[1]) == GenSym 
        symbolInfo[ast.args[1]] = typeof(ast.args[2])
        dprintln(2,"Add symbolInfo for GenSym: ", ast.args[1], "==>", typeof(ast.args[2]))
    end
    return nothing
end

function initSymbol2TypeDict(expr)
    assert(expr.head == :lambda) # (:lambda, {param, meta@{localvars, types, freevars}, body})

    local varinfo = expr.args[2][2]
    symbolInfo = Dict{Union(Symbol, Integer), Any}()
    for i = 1:length(varinfo)
        symbolInfo[varinfo[i][1]] = varinfo[i][2]
        dprintln(4,"Add symbolInfo: ", varinfo[i][1], "==>", symbolInfo[varinfo[i][1]])
    end
    
    # Add GenSym's types. We do not have to walk the AST. Instead, lambda has that info
    # CompilerTools.AstWalker.AstWalk(expr.args[3], memoize_GenSym_types, symbolInfo)
    local gensym_types = expr.args[2][4]
    for id = 1:length(gensym_types)
        symbolInfo[id - 1] = gensym_types[id] # GenSym id starts from 0
        dprintln(4,"Add GenSym symbolInfo: ", id - 1, "==>", symbolInfo[id - 1])
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
  elseif asttyp == GenSym
    return get(symbolInfo, ast.id, Nothing)
  else
    return asttyp
  end
end

function analyze_type(typ)
    is_number = (typ <: Number)
    # ASSUMPTION: we assume the user program does not use range 
    # for any array computation, although Range is a subtype of AbstractArray
    is_array  = (typ <: AbstractArray && !(typ <: Range))
    is_number, is_array
end

function analyze_types(result_type, arg_types :: Tuple)
    all_numbers, some_arrays = analyze_type(result_type)
    for t in arg_types
        is_number, is_array = analyze_type(t)
        all_numbers =  all_numbers && is_number
        some_arrays = some_arrays || is_array
    end
    all_numbers, some_arrays
end

function name_of_module_or_function(arg)
    if typeof(arg) == Symbol
        return string(arg)
    elseif typeof(arg) == QuoteNode && typeof(arg.value) == Symbol
        return string(arg.value)
    else
        throw(UnhandledModuleOrFunctionName(arg))
    end
end

function resolve_module_function_names(call_args)
    assert(length(call_args) > 0)
    module_name, function_name = "", ""
    if typeof(call_args[1]) == Symbol # Example: :*
        function_name = string(call_args[1])
    elseif isa(call_args[1], TopNode) && # Example: top(getfield), SparseAccelerator,:SpMV
            length(call_args) == 3 
        if call_args[1] == TopNode(:getfield)
            module_name = name_of_module_or_function(call_args[2])
            function_name = name_of_module_or_function(call_args[3])
        else
            function_name = name_of_module_or_function(call_args[1].name)
        end
    elseif  isa(call_args[1], Expr) &&
            call_args[1].head == :call # Example: (:call, top(getfield), SparseAccelerator,:SpMV)
        return resolve_module_function_names(call_args[1].args)
    end
    module_name, function_name
end

function checkDistributivityForCall(expr, symbolInfo, distributive)
    head = expr.head
    args = expr.args
        
    # A typical call are in the following forms
    #   Expr(:call, :*, :A, :x)
    #   Expr(:call, :(:call, top(getfield), SparseAccelerator,:SpMV), :A, :x)
    # So the first argument is the function, the others are the arguments for it.
        
    arg_types = ntuple(i-> typeOfNode(args[i+1], symbolInfo), length(args) - 1)
    all_numbers, some_arrays = analyze_types(expr.typ, arg_types)
    if all_numbers || !some_arrays
        # Result and args are all numbers, or there may be other types (like 
        # Range{UInt64}) but no regular arrays. This function is distributive.
        # Do nothing
    else
        module_name, function_name = resolve_module_function_names(args)
        if function_name == ""
            throw(UnresolvedFunction(head, args[1]))
        end
        fd = lookup_function_description(module_name, function_name, arg_types)
        if fd != nothing
            distributive = distributive && fd.distributive
        else
            throw(UndescribedFunction(module_name, function_name, arg_types))
        end
    end
    
    if distributive
        # In case some args are trees themselves (even if the args' result types are
        # Number), we should check the args further
        for x in args[2 : end]
            distributive = checkDistributivity(x, symbolInfo, distributive)
        end
        # ISSUE: what if an arg is an array element? It is a scalar, but it involves an array.
    end
    return distributive
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
function checkDistributivity(ast::Any, symbolInfo::Dict, distributive::Bool)
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
    elseif is_assignment(ast)  # Updating operators
        assert(length(args) == 2)
        local lhs = args[1]
        local rhs = args[2]
        distributive = checkDistributivity(rhs, symbolInfo, distributive)
        dprintln(2, "\tAssignment ", distributive)
    elseif head == :return
        distributive = checkDistributivity(args, symbolInfo, distributive)
        dprintln(2, "\tReturn ", distributive)
    elseif head == :call || head == :call1
        distributive = checkDistributivityForCall(ast, symbolInfo, distributive)
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
  elseif asttyp == Symbol || asttyp == GenSym
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

      # There can be such arrays like [:(((top(getfield))(SparseAccelerator,:SpMV))(A,x)]
      # So we have to check each element of the array.
      for element in  ast
          distributive = checkDistributivity(element, symbolInfo, distributive)
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
  loop_info = CompilerTools.Loops.compute_dom_loops(lives.cfg)

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
function SpMV!(w::Vector, alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector, gamma::Number)
    assert(length(w) == length(y))

    if DEFAULT_LIBRARY == PCL_LIB
        A2 = CreateCSR(A)
        ccall((:CSR_MultiplyWithVector, LIB_PATH), Void,
              (Ptr{Cdouble}, Cdouble, Ptr{Void}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cdouble),
              pointer(w), alpha, A2, pointer(x), beta, pointer(y), gamma)
        DestroyCSR(A2)
    else
        # use Julia implementation
        w = alpha * A * x + beta * y + gamma
    end
end

# y = A*x
SpMV!(y::Vector, A::SparseMatrixCSC, x::Vector) = SpMV!(y, one(eltype(x)), A, x, zero(eltype(y)), y, zero(eltype(x)))

function PageRank!(w::Vector, alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector, gamma::Number, z::Vector)
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
function SpMV(alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector, gamma::Number)
  w = Array(Cdouble, length(x))
  SpMV!(w, alpha, A, x, beta, y, gamma)
  w
end

# alpha*A*x + beta*y
SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector, beta::Number, y::AbstractVector) = SpMV(alpha, A, x, beta, y, zero(eltype(x)))

# alpha*A*x + y
SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector, y::AbstractVector) = SpMV(alpha, A, x, one(eltype(y)), y, zero(eltype(x)))

# alpha*A*x
SpMV(alpha::Number, A::SparseMatrixCSC, x::AbstractVector) = SpMV(alpha, A, x, zero(eltype(x)), x, zero(eltype(x)))

# A*x
SpMV(A::SparseMatrixCSC, x::Vector) = SpMV(one(eltype(x)), A, x, zero(eltype(x)), x, zero(eltype(x)))

function SpMV_conditional_reordering(A::SparseMatrixCSC, x::Vector, reorder_done, beneficial)
    if reorder_done
        return SpMV(A, x)
    end
    time1 = time()
    y = SpMV(A, x)
    time2 = time()
    nnz = size(A.nzval, 1)
    rows = A.m
    SpMV_bandwidth = (nnz * 12 + rows * 3 * 8) / (time2 - time1) * 128 / 10 ^ 9
    machine_peak_bandwidth = 60
    if abs(SpMV_bandwidth - machine_peak_bandwidth) > 15
        beneficial[1] = true
    end
    return y
end

function init_conditional_reordering(beneficial, reorder_done)
    # Benefical can be just a scalar var. But we will treat it as a 1-element array
    # so that we do not have difficulty in changing it in calling 
    # SpMV_conditional_reordering
    beneficial = Vector{Bool}()
    push!(beneficial, false)
    reorder_done = false
end

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
  w
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
  w
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
  w
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

#function Base.A_mul_B!(alpha::Number, A::SparseMatrixCSC, x::Vector, beta::Number, y::Vector)
#  SparseAccelerator.SpMV(alpha, A, x, beta, y)
#end

