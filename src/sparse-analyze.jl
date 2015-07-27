export CSR_ReorderMatrix, reorderVector, reverseReorderVector

# This controls the debug print level.  0 prints nothing.  At the moment, 2 prints everything.
#DEBUG_LVL=0

#function set_debug_level(x)
#    global DEBUG_LVL = x
#end

# A debug print routine.
#function dprint(level,msgs...)
#    if(DEBUG_LVL >= level)
#        print(msgs...)
#    end 
#end

# A debug print routine.
#function dprintln(level,msgs...)
#    if(DEBUG_LVL >= level)
#        println(msgs...)
#    end 
#end

const LIB_PATH = libcsr

function fwdTriSolve!(A::SparseMatrixCSC, B::AbstractVecOrMat, fknob)
  ccall((:ForwardTriangularSolve, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(B), pointer(B), fknob == nothing ? C_NULL : pointer(fknob))
  return B
end
    
#function new_matrix_knob
#  return ccall((:NewMatrixKnob, LIB_PATH), Ptr{Void}, ())
#end

#function increment_matrix_version
#  return ccall((:IncrementMatrixVersion, LIB_PATH), Ptr{Void}, ())
#end

function delete_matrix_knob(mknob)
    ccall((:DeleteMatrixKnob, LIB_PATH), Void, (Ptr{Void},), mknob)
end

function add_matrix_knob(fknob, mknob)
  ccall((:AddMatrixKnob, LIB_PATH), Void, (Ptr{Void}, Ptr{Void}), 
    pointer(fknob), pointer(mkob))
end

function new_function_knob(fknob_creator)
  return ccall(fknob_creator, Ptr{Void}, ())
end

function delete_function_knob(fknob_deletor, fknob)
    ccall(fknob_deletor, Void, (Ptr{Void},), fknob)
end

# In reordering, we insert some calls to the following 3 functions. So they are executed secretly
# Reorder sparse matrix A and store the result in newA. A itself is not changed.
function CSR_ReorderMatrix(A :: SparseMatrixCSC, 
                           newA :: SparseMatrixCSC, 
                           P :: Vector, 
                           Pprime :: Vector, 
                           getPermutation :: Bool, 
                           oneBasedInput :: Bool, 
                           oneBasedOutput :: Bool)
  ccall((:CSR_ReorderMatrix, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool, Bool, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(newA.colptr), pointer(newA.rowval), pointer(newA.nzval),
               pointer(P), pointer(Pprime), getPermutation, oneBasedInput, oneBasedOutput)
end

function CSR_Bandwidth(A::SparseMatrixCSC)
   A2 = CreateCSR(A)
   bw = ccall((:CSR_GetBandwidth, LIB_PATH), Cint,
         (Ptr{Void},),
         A2)
   DestroyCSR(A2)
   bw
end

function reorderVector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVector, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(newV), pointer(P), length(V))
end

function reverseReorderVector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVectorWithInversePerm, LIB_PATH), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(newV), pointer(P), length(V))
end

function allocateForPermutation(M::Symbol, new_stmts)
    P = gensym("P")
    stmt = :($P = Array(Cint, size($M, 2)))
    push!(new_stmts, stmt)
    
    Pprime = gensym("Pprime")
    stmt = :($Pprime = Array(Cint, size($M, 2)))
    push!(new_stmts, stmt)
    (P, Pprime)
end

# A is the matrix to be reordered. If getPermutation is true, we also compute 
# permutation and inverse permutation
# info (P and Pprime). Otherwise, the info has already been computed. That means, this
# function must be called to reorder matrix M first, and then for other matrices
function reorderMatrix(A::Symbol, P::Symbol, Pprime::Symbol, new_stmts, getPermuation::Bool, oneBasedInput::Bool, oneBasedOutput::Bool)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the sparse matrix format and element type. Should make it
    # general in future
    newA = gensym(string(A))
    stmt = Expr(:(=), newA,
        Expr(:call, :SparseMatrixCSC, 
                  Expr(:call, TopNode(:getfield), A, QuoteNode(:m)),
                  Expr(:call, TopNode(:getfield), A, QuoteNode(:n)), 
                  Expr(:call, :Array, :Cint, 
                    Expr(:call, TopNode(:arraylen), 
                        Expr(:call, TopNode(:getfield), A, QuoteNode(:colptr)))),
                  Expr(:call, :Array, :Cint, 
                    Expr(:call, TopNode(:arraylen), 
                        Expr(:call, TopNode(:getfield), A, QuoteNode(:rowval)))),
                  Expr(:call, :Array, :Cdouble, 
                    Expr(:call, TopNode(:arraylen), 
                        Expr(:call, TopNode(:getfield), A, QuoteNode(:nzval))))                        
        )
    )
    push!(new_stmts, stmt)
    
    # Do the actual reordering in the C library
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :CSR_ReorderMatrix), A, newA, P, Pprime, getPermuation, oneBasedInput, oneBasedOutput)
    push!(new_stmts, stmt)
    
    # Update the original matrix with the new data. Note: assignment between
    # two arrays actually make them aliased. So there is no copy cost
    stmt = Expr(:(=), A, newA)
    push!(new_stmts, stmt)
end

reverseReorderMatrix(sym, P, Pprime, landingPad, getPermuation, oneBasedInput, oneBasedOutput) = 
       reorderMatrix(sym, Pprime, P, landingPad, getPermuation, oneBasedInput, oneBasedOutput)

function reorderVector(V::Symbol, P::Symbol, new_stmts)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the element type. Should make it general in future
    newV = gensym(string(V))
    stmt = Expr(:(=), newV,
                Expr(:call, :Array, :Cdouble, 
                    Expr(:call, TopNode(:arraylen), V)))
    push!(new_stmts, stmt)
    
    # Do the actual reordering in the C library
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :reorderVector),
                V, newV, P)
    push!(new_stmts, stmt)
    
    # Update the original vector with the new data
    stmt = Expr(:(=), V, newV )
    push!(new_stmts, stmt)
end

function reverseReorderVector(V::Symbol, P::Symbol, new_stmts)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the element type. Should make it general in future
    newV = gensym(string(V))
    stmt = Expr(:(=), newV,
                Expr(:call, :Array, :Cdouble, 
                    Expr(:call, TopNode(:arraylen), V)))
    push!(new_stmts, stmt)
    
    # Do the actual reordering in the C library
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :reverseReorderVector),
                V, newV, P)
    push!(new_stmts, stmt)
    
    # Update the original vector with the new data
    stmt = Expr(:(=), V, newV )
    push!(new_stmts, stmt)
end

function insertTimerBeforeLoop(new_stmts)
    t1 = gensym("t1")
    stmt = Expr(:(=), t1,
                Expr(:call, TopNode(:ccall), QuoteNode(:clock_now), :Float64, 
                    Expr(:call1, TopNode(:tuple))
                )
            )
    push!(new_stmts, stmt)
    t1
end

function insertTimerAfterLoop(t1, new_stmts)
    t2 = gensym("t2")
    stmt = Expr(:(=), t2,
                Expr(:call, TopNode(:ccall), QuoteNode(:clock_now), :Float64, 
                    Expr(:call1, TopNode(:tuple))
                )
            )
    push!(new_stmts, stmt)

    t3 = gensym("t3")
    stmt = Expr(:(=), t3,
                Expr(:call, TopNode(:box), :Float64, 
                    Expr(:call, TopNode(:sub_float), t2, t1)
                )
            )
    push!(new_stmts, stmt)

    stmt = Expr(:call, :println, GlobalRef(Base, :STDOUT), "Time of loop= ", t3, " seconds")
    push!(new_stmts, stmt)
end

# Old version, not differentiating LHS and RHS symbols. To remove
function addToIA(currentNode, IA, symbolInfo)
    if typeOfNode(currentNode, symbolInfo) <: AbstractArray
        if typeof(currentNode) == Symbol
            push!(IA, currentNode) 
        elseif typeof(currentNode) == SymbolNode
            push!(IA, currentNode.name)
        end
    end 
end

function printSpaces(i)
    for j = 1:i
        print("\t")
    end
end

# Old version, not differentiating LHS and RHS symbols. To remove
function findIA(currentNode, IA, seeds, symbolInfo, i)
#   printSpaces(i)
#   println("findIA: ", currentNode)
    addToIA(currentNode, IA, symbolInfo)
    if typeof(currentNode) <: Expr
        for c in currentNode.args
#            printSpaces(i+1)
#            println("c=", c, " type=", typeOfNode(c, symbolInfo))
            if typeOfNode(c, symbolInfo) <: Number
                push!(seeds, c)
            else
                findIA(c, IA, seeds, symbolInfo, i+1)
            end
        end
    end
end

function deepType(x)
    if typeof(x) == SymbolNode
        return x.typ
    elseif typeof(x) == Expr
        return x.typ
    else
        return typeof(x)
    end
end

function isMul(x)
  if typeof(x) == Expr && x.head == :call && x.args[1] == :(*)
    lhs = x.args[2]
    rhs = x.args[3]
    if deepType(lhs) <: Number && deepType(rhs) <: Vector
      return (true,lhs,rhs)
    elseif deepType(rhs) <: Number && deepType(lhs) <: Vector
      return (true,rhs,lhs)
    end
  end
  return (false,nothing,nothing)
end

function isWaxpby(node)
  if typeof(node)         == Expr && node.head         == :call && 
     typeof(node.args[1]) == Expr && node.args[1].head == :call && 
     node.args[1].args[1] == TopNode(:getfield) && node.args[1].args[2] == :SparseAccelerator && node.args[1].args[3] == QuoteNode(:WAXPBY)
    return (true, getSName(node.args[3]), getSName(node.args[5]))
  end
  return (false, nothing, nothing)
end

function isSpMV(node)
  if typeof(node)         == Expr && node.head         == :call && 
     typeof(node.args[1]) == Expr && node.args[1].head == :call && 
     node.args[1].args[1] == TopNode(:getfield) && node.args[1].args[2] == :SparseAccelerator && node.args[1].args[3] == QuoteNode(:SpMV)
    return true
  end 
  return false
end

function getSName(ssn)
  stype = typeof(ssn)
  if stype == Symbol
    return ssn
  elseif stype == SymbolNode
    return ssn.name
  elseif stype == Expr && ssn.head == :(::)
    return ssn.args[1]
  end

  dprintln(0, "getSName ssn = ", ssn, " stype = ", stype)
  if stype == Expr
    dprintln(0, "ssn.head = ", ssn.head)
  end
  throw(string("getSName called with something of type ", stype))
end

function optimize_calls(ast, state, top_level_number, is_top_level, read)
  asttyp = typeof(ast)
  dprintln(3,"optimize_calls ast = ", ast, " type = ", typeof(ast))
  if asttyp == Expr && ast.head == :call
    dprintln(3,"optimize_calls found call expr ", ast)
    if ast.args[1] == :A_mul_B!
      dprintln(3,"optimize_calls found A_mul_B!")
      if length(ast.args) == 4
        y = ast.args[2]
        A = ast.args[3]
        x = ast.args[4]
        if A.typ <: SparseMatrixCSC && x.typ <: Vector
          ast.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!))
          return ast
        end
      elseif length(ast.args) == 6
        alpha = ast.args[2]
        A = ast.args[3]
        x = ast.args[4]
        beta = ast.args[5]
        y = ast.args[6]
        if A.typ <: SparseMatrixCSC && x.typ <: Vector
          ast.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!))
          return ast
        end
      end
    elseif ast.args[1] == :dot
      dprintln(3,"optimize_calls found :dot ", ast)
      assert(length(ast.args) == 3)
      ast.args[2] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[2], optimize_calls, nothing))
      ast.args[3] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[3], optimize_calls, nothing))
      arg1 = ast.args[2]
      arg2 = ast.args[3]
      dprintln(3,"arg1 = ", arg1, " arg2 = ", arg2, " arg1.typ = ", deepType(arg1), " arg2.typ = ", deepType(arg2))
      if deepType(arg1) <: Vector && deepType(arg2) <: Vector
        dprintln(3,"optimize_calls converting to SparseAccelerator.Dot")
        ast.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:Dot))
      end
      return ast
    elseif ast.args[1] == :(+) || ast.args[1] == :(-)
      if length(ast.args) > 2   # exclude processing for unary -
        dprintln(3,"optimize_calls found + or -")
        for i = 2:length(ast.args)
          ast.args[i] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[i], optimize_calls, nothing))
        end
        num_operands = length(ast.args) - 1
        num_replaced = 0
        for rawi = 2 : num_operands
          i = rawi - num_replaced
          arg1 = ast.args[i]
          arg2 = ast.args[i+1]
          dprintln(3,"arg1 = ", arg1, " arg2 = ", arg2, " arg1.typ = ", deepType(arg1), " arg2.typ = ", deepType(arg2))
          if deepType(arg1) <: Vector && deepType(arg2) <: Vector 
            (is_mul1, scalar1, vec1) = isMul(arg1)
            (is_mul2, scalar2, vec2) = isMul(arg2)
            dprintln(3,"optimize_calls found vector/vector +. mul1 = ", is_mul1, " mul2 = ", is_mul2)
            if is_mul1 || is_mul2
              replacement = Expr(:call)
              replacement.args = Array(Any,5)
              dprintln(3,"optimize_calls converting to SparseAccelerator.WAXPBY")
              replacement.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY))
              if is_mul1
                replacement.args[2] = scalar1
                replacement.args[3] = vec1
              else
                replacement.args[2] = 1
                replacement.args[3] = arg1
              end
              if is_mul2
                if ast.args[1] == :(+)
                  replacement.args[4] = scalar2
                else
                  replacement.args[4] = CompilerTools.LivenessAnalysis.TypedExpr(deepType(scalar2), :call, :(-), scalar2)
                end
                replacement.args[5] = vec2
              else
                if ast.args[1] == :(+)
                  replacement.args[4] = 1
                else
                  replacement.args[4] = -1
                end
                replacement.args[5] = arg2
              end
              ast.args = [ast.args[1:i-1]; Expr(:call, CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)), arg1, arg2); ast.args[i+2:length(ast.args)]]
              num_replaced = num_replaced + 1
            elseif isSpMV(arg1) && length(arg1.args) == 4 # 4 = SpMV of the form alpha*A*x
              if ast.args[1] == :(+)
                ast.args = [ast.args[1:i-1]; Expr(:call, CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)), arg1.args[2], arg1.args[3], arg1.args[4], arg2); ast.args[i+2:length(ast.args)]]
              else
                # Do I need to type "-1" in some way here?
                ast.args = [ast.args[1:i-1]; Expr(:call, CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)), arg1.args[2], arg1.args[3], arg1.args[4], -1, arg2); ast.args[i+2:length(ast.args)]]
              end
              num_replaced = num_replaced + 1
            end
          end
        end
        if length(ast.args) == 2
          dprintln(3,"optimize_calls up-leveling the SpMV call")
          ast.args = ast.args[2].args
        end
        return ast
      end
    elseif ast.args[1] == :(*)
      dprintln(3,"optimize_calls found *")
      for i = 2:length(ast.args)
        ast.args[i] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[i], optimize_calls, nothing))
      end
      num_operands = length(ast.args) - 1
      for i = num_operands : -1 : 2
        arg1 = ast.args[i]
        arg2 = ast.args[i+1]
        dprintln(3,"arg1 = ", arg1, " arg2 = ", arg2, " arg1.typ = ", deepType(arg1), " arg2.typ = ", deepType(arg2))
        if deepType(arg1) <: SparseMatrixCSC && deepType(arg2) <: Vector
          ast.args = [ast.args[1:i-1]; Expr(:call, CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)), arg1, arg2); ast.args[i+2:length(ast.args)]]
          dprintln(3,"optimize_calls converting to 2 argument SparseAccelerator.SpMV, length = ", length(ast.args))
        elseif deepType(arg1) <: Number && isSpMV(arg2) && length(arg2.args) == 3  # 3 = SpMV + two arguments
          ast.args = [ast.args[1:i-1]; Expr(:call, CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)), arg1, arg2.args[2], arg2.args[3]); ast.args[i+2:length(ast.args)]]
          dprintln(3,"optimize_calls converting to 3 argument SparseAccelerator.SpMV, length = ", length(ast.args))
        end
      end 
      # If the length of ast.args has been shortened to 2 then we can no longer use the "*" node so we elevate ast.args[1] up one level.
      if length(ast.args) == 2
        dprintln(3,"optimize_calls up-leveling the SpMV call")
        ast.args = ast.args[2].args
      end
      return ast
    end
  elseif asttyp == Expr && ast.head == :(=)
    dprintln(3,"optimize_calls found = ", ast)
    ast.args[1] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[1], optimize_calls, nothing))
    ast.args[2] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[2], optimize_calls, nothing))
    dprintln(3,"after recursive optimization ", ast)
    (rhs_waxpby, rhs_x, rhs_y) = isWaxpby(ast.args[2])
    dprintln(3, "rhs_waxpby = ", rhs_waxpby, " ", rhs_x, " ", rhs_y)
    if rhs_waxpby && (ast.args[1] == rhs_x || ast.args[1] == rhs_y)
      dprintln(3,"optimize_calls converting to SparseAccelerator.WAXPBY!")
      ast.args[2].args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY!))
      splice!(ast.args[2].args, 2:5, [ast.args[1]; ast.args[2].args[2:5]])
      ast = ast.args[2]
    end
    return ast
  end
  return nothing
end

# Try to reorder sparse matrices and related (dense vectors) in the given loop
# L is a loop region. M is a sparse matrix live into L. 
# Lives is the set of matrices or vectors live into/out of L
# ASSUMPTION: if two matrices/vectors alias, they must be fully overlap.
#             We cannot handle partial overlapping aliases.
# ISSUE: even for this simple full-or-none-alias assumption, we have
#        difficult to analyze aliases at compile time. And whether two function
#        parameters alias or not is dependent on the specific call site
#        at run time. So we would simply assume NO matrix/vector alias 
#        with any other one.
# TODO:  check that no two matrices/vectors alias dynamically at runtime
#        by inserting the check to the function AST. The code would be like
#        this: 
#            if (A==B) aliased=true
#            if (!aliased) M=P'*M*P
# TODO:  consider matrix sizes and types. We need to make sure that all matrices
#        have the same size and type(?)
function reorderLoop(L, M, lives, symbolInfo)
    dprintln(2, "Reorder: loop=", L, " Matrix=", M)

    if(DEBUG_LVL >= 2)
        println("******** CFG before reordering: ********")
        show(lives.cfg);
    end 

    # If we reorder M inside L, consequently, some other matrices or vectors
    # need to be reordered. For example, for the following code
    #       N = M * v 
    #       K = N + speye(size(N, 1))
    # v has to be reordered too, to make M * v meaningful. As the result of 
    # M * v, N got reordered as well. Consequently, the matrix due to speye()
    # has to be reordered, and K as well. This is a propagation.
    # Assumption: we know the semantics of all the operators/functions.

    bbs = lives.basic_blocks
    
    # Build inter-dependent arrays
    IAs = Dict{Any, Set}()    
    for bbnum in L.members
        for stmt_index = 1:length(bbs[lives.cfg.basic_blocks[bbnum]].statements)
            # Replace calls to optimized versions provided in SparseAccelerator module.
            bbs[lives.cfg.basic_blocks[bbnum]].statements[stmt_index].expr = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(bbs[lives.cfg.basic_blocks[bbnum]].statements[stmt_index].tls.expr, optimize_calls, nothing))
            stmt = bbs[lives.cfg.basic_blocks[bbnum]].statements[stmt_index]
            IAs[stmt] = Set{Any}()
            seeds = Set{Any}()
            push!(seeds, stmt.expr)
            while !isempty(seeds)
                seed = pop!(seeds)
                IA = Set{Any}()
                findIA(seed, IA, seeds, symbolInfo, 1)
                if !isempty(IA)
                    push!(IAs[stmt], IA)
                end
            end
#            println("###Stmt: ", stmt.expr)
#            println("  IAs:", IAs[stmt])        
        end
    end    

    reordered = Set{Symbol}()
    push!(reordered, M)

    changed = true
    while changed
        changed = false
        for bbnum in L.members
            for stmt in bbs[lives.cfg.basic_blocks[bbnum]].statements
                for IA in IAs[stmt]
                    if ⊈(IA, reordered)
                        if !isempty(intersect(IA, reordered))
                            union!(reordered, IA)
                            changed = true
                        end
                    end
                end
            end
        end
    end
    
    dprintln(2, "Reordered:", reordered)
 
    # The symbols that should be reordered before L are the reorderedUses live into L
    headBlock = bbs[L.head]
    reorderedBeforeL = intersect(reordered, headBlock.live_in)
    
    dprintln(2, "To be reordered before L: ", reorderedBeforeL)

    if isempty(reorderedBeforeL)
        return
    end
    
    #TODO: benefit-cost analysis

    # Only those arrays that are updated in the loop need to be reverse reordered afterwards
    updatedInLoop = Set{Symbol}()
    for bbnum in L.members
        union!(updatedInLoop, bbs[bbnum].def)
    end
    dprintln(2, "updatedInLoop: ", updatedInLoop)
      
    # New a vector to hold the new reordering statements R(LiveIn) before the loop. 
    # We do not insert them into the CFG at this moment, since that will change the
    # pred-succ and live info, leading to some subtle errors. We need CFG not changed
    # until all new statements are ready.
    new_stmts_before_L = Expr[]

    # Allocate space to store the permutation and inverse permutation info
    (P, Pprime) = allocateForPermutation(M, new_stmts_before_L)

    # Compute P and Pprime, and reorder M
    reorderMatrix(M, P, Pprime, new_stmts_before_L, true, true, true)
    
    # Now reorder other arrays
    for sym in reorderedBeforeL
        if sym != M
            if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                reorderMatrix(sym, P, Pprime, new_stmts_before_L, false, true, true)
            else
                reorderVector(sym, P, new_stmts_before_L)
            end
        end
    end
    
    # For perf measurement only. 
    # TODO: add a switch here
#    t1 = insertTimerBeforeLoop(new_stmts_before_L)
             
    # At the exit of the loop, we need to reverse reorder those that live out of the loop,
    # to recover their original order before they are getting used outside of the loop
    # We remember those in an array. Each element is a tuple 
    # (bb label, succ label, the new statements to be inserted on the edge from bb to succ)
    new_stmts_after_L =  Any[]
    reorderedAndUpdated = intersect(reordered, updatedInLoop)
    for bbnum in L.members
        bb = bbs[bbnum]
        for succ in bb.succs
            if !in(succ.label, L.members)
                # For perf measurement only
                # TODO: add a switch here
                new_stmts = (bbnum, succ.label,  Expr[])
                push!(new_stmts_after_L, new_stmts)
#                insertTimerAfterLoop(t1, new_stmts[3])

                reverseReordered = intersect(reorderedAndUpdated, succ.live_in)
                if isempty(reverseReordered)
                    continue
                end
                dprintln(2, "ReverseReorder on edge ", bbnum, " --> ", succ.label)
                
                for sym in reverseReordered
                    if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                        reverseReorderMatrix(sym, P, Pprime, new_stmts[3], false, true, true)
                    else
                        reverseReorderVector(sym, P, new_stmts[3])
                    end
                end
            end
        end
    end

    # Now actually change the CFG.
    (new_bb, new_goto_stmt) = CompilerTools.CFGs.insertBefore(lives.cfg, L.head, true, L.back_edge)
    lives.basic_blocks[new_bb] = CompilerTools.LivenessAnalysis.BasicBlock(new_bb)
    for stmt in new_stmts_before_L
        CompilerTools.CFGs.addStatementToEndOfBlock(lives.cfg, new_bb, stmt)
        # Do we need to create the statement in the LivenessAnalysis block as well?
    end
    if new_goto_stmt != nothing
      push!(new_bb.statements, new_goto_stmt)
      # Do we need to create the statement in the LivenessAnalysis block as well?
    end
    
    for (pred, succ, new_stmts) in new_stmts_after_L
        (new_bb, new_goto_stmt) = CompilerTools.CFGs.insertBetween(lives.cfg, pred, succ)
        lives.basic_blocks[new_bb] = CompilerTools.LivenessAnalysis.BasicBlock(new_bb)
        for stmt in new_stmts
            CompilerTools.CFGs.addStatementToEndOfBlock(lives.cfg, new_bb, stmt)
        end
        if new_goto_stmt != nothing
          push!(new_bb.statements, new_goto_stmt)
        end
    end

    if(DEBUG_LVL >= 2)
        println("******** CFG after reordering: ********")
        show(lives.cfg);
    end 
end

# Map each BB to a Bool: true if the BB in any loop 
function BBsInLoop(lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                   loop_info :: CompilerTools.Loops.DomLoops)
    in_loop = Dict{Int, Bool}()
    for (j,bb) in lives.cfg.basic_blocks
        in_loop[bb.label] = false
    end

    for L in loop_info.loops
        for bbnum in L.members
            in_loop[bbnum] = true
        end
    end
    return in_loop
end

# Map each statement to a Bool: true if the statement in any loop 
function statementsInLoop(lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                   loop_info :: CompilerTools.Loops.DomLoops)
    BB_in_loop = BBsInLoop(lives, loop_info)
    stmt_BB = Dict{CompilerTools.CFGs.TopLevelStatement, Int}()
    for (j, bb) in lives.cfg.basic_blocks
        for stmt in bb.statements
            stmt_BB[stmt] = bb.label
        end
    end
    return BB_in_loop, stmt_BB
end

function findMmread(lives :: CompilerTools.LivenessAnalysis.BlockLiveness, in_loop :: Dict{Int,Bool})
    # Look for such a statement outside any loop: e.g. A = ((top(getfield))(MatrixMarket,:mmread))("b.mtx")
    for (j,bb) in lives.cfg.basic_blocks
        if in_loop[bb.label] 
            continue
        end
        for stmt_idx = 1 : length(bb.statements)
            expr = bb.statements[stmt_idx].expr
            if typeof(expr) != Expr
                continue
            end
            if expr.head == :(=) 
                lhs  = expr.args[1]
                if typeof(lhs) != Symbol && typeof(lhs) != GenSym
                    continue
                end
                rhs  = expr.args[2]
                if typeof(rhs) == Expr && rhs.head == :call && size(rhs.args, 1) == 2 #two arguments: ((top(getfield))(MatrixMarket,:mmread)), and filename 
                    func = rhs.args[1] #:((top(getfield))(MatrixMarket,:mmread))
                    if typeof(func) == Expr && func.head == :call && size(func.args, 1) == 3
                        if func.args[1] == TopNode(:getfield) &&  
                            func.args[2] == :MatrixMarket && 
                            func.args[3] == QuoteNode(:mmread) # Note: some difference between interactive env and comand line. In interactive env, QuoteNode(:mmread) should be written as ":mmread" to take effect                            
                            return bb, stmt_idx
                        end
                    end
                end
            end
        end
    end
    return nothing, 0
end 

type ExitEdge
    from_BB       :: CompilerTools.LivenessAnalysis.BasicBlock
    from_stmt_idx :: Int
    to_BB         :: CompilerTools.LivenessAnalysis.BasicBlock
    to_stmt_idx   :: Int
end

# One part of a region. Both from/to statements and the other statements between them are included
type RegionInterval
    BB            :: CompilerTools.LivenessAnalysis.BasicBlock
    from_stmt_idx :: Int
    to_stmt_idx   :: Int
    succs         :: Set{RegionInterval}
    preds         :: Set{RegionInterval}
    
    RegionInterval(bb, from, to) = new(bb, from, to, Set{RegionInterval}(), Set{RegionInterval}())
end

type Region
    first_BB        :: CompilerTools.LivenessAnalysis.BasicBlock # the first BB that dominates all other BBs in the region
    mmread_stmt_idx :: Int # 0 unless first_BB has a statement containing mmread(), in which case this points to the index of the statement
    exits           :: Set{ExitEdge}
    intervals       :: Vector{RegionInterval}
end

module Base
    function show_expr_type(io::IO, ty)
        return
    end
end

function show_interval(interval :: RegionInterval, level)
    for i = 1:level print("  ") end
    println("BB ", interval.BB.cfgbb.label, " stmts ", interval.from_stmt_idx, "(", 
        1 <= interval.from_stmt_idx && interval.from_stmt_idx <= length(interval.BB.statements) ?
        interval.BB.statements[interval.from_stmt_idx].tls.index : "", ") : ", interval.to_stmt_idx, "(",
        1 <= interval.to_stmt_idx && interval.to_stmt_idx <= length(interval.BB.statements) ?
        interval.BB.statements[interval.to_stmt_idx].tls.index : "", ")")
end

# Note: to be sure of profitability, we should require each path to pass through a loop. 
# However, in pagerank-mmread-inside.jl pagerank(), the CFG is like this:
# Block -1: Succ( -3 3 )    
#           mmread 
#           repeat = 100
#           GenSym(0) = colon(1,repeat::Int64)::UnitRange{Int64}
#           s135 = (top(start))(GenSym(0))::Int64
#           unless (top(!))((top(done))(GenSym(0),#s135::Int64)::Bool)::Bool goto 3
# Block -3: has a loop, and finally toward Block 3
# Block  3: return 
# So there are totally two paths. One of them (-1-->3) does not have any loop.
# This is because we do not have analysis to tell us that -1-->3 is an infeasible path.
# For now, give up with the requirement of having a loop.
# TODO: enforce the requirement after some analysis like DCE is implemented
function DFSGrowRegion(first_BB :: CompilerTools.LivenessAnalysis.BasicBlock, 
                       current_BB :: CompilerTools.LivenessAnalysis.BasicBlock, 
                       start_stmt_idx :: Int64, 
                       lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                       in_loop :: Dict{Int,Bool}, 
                       visited :: Dict{Int,Bool}, 
                       has_loop :: Dict{Int,Bool}, 
                       bb_interval :: Dict{Int,RegionInterval}, 
                       exits :: Set{ExitEdge}, 
                       intervals :: Vector{RegionInterval}, 
                       symbolInfo :: Dict{Union(Symbol,Integer),Any}, 
                       loop_info :: CompilerTools.Loops.DomLoops, 
                       loop_bbs)
    if (DEBUG_LVL >= 2)
        println("\n\nDFSGrowRegion from BB ", current_BB.cfgbb.label, " stmt ", start_stmt_idx, "(", 
            1 <= start_stmt_idx && start_stmt_idx <= length(current_BB.statements) ? 
            current_BB.statements[start_stmt_idx].tls.index : "", ")");
    end
    
    # Disable the requirement of having a loop.
    # TODO: enable it in future
    loop_expected_on_every_path = false
    
    assert(!visited[current_BB.cfgbb.label])
    visited[current_BB.cfgbb.label]  = true
    has_loop[current_BB.cfgbb.label] = in_loop[current_BB.cfgbb.label]
    last_stmt_idx = length(current_BB.cfgbb.statements)
    for stmt_idx = start_stmt_idx : last_stmt_idx
        expr = current_BB.cfgbb.statements[stmt_idx].expr
        if typeof(expr) != Expr
            continue
        end
        if expr.head == :return
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, nothing
            end
            exit = ExitEdge(current_BB, stmt_idx - 1, current_BB, stmt_idx)
            push!(exits, exit)
            interval = RegionInterval(current_BB, start_stmt_idx, stmt_idx - 1)
            push!(intervals, interval)
            bb_interval[current_BB.cfgbb.label] = interval
            return true, interval
        end
        if expr.head == :throw
            throw("DFSGrowRegion: throw not handled")
            return false, nothing
        end
        distributive = checkDistributivity(expr, symbolInfo, true)
        if !distributive
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, nothing
            end
            if loop_bbs != nothing
                # The region cannot cover the whole loop
                return false, nothing
            end
            exit = ExitEdge(current_BB, stmt_idx - 1, current_BB, stmt_idx)
            push!(exits, exit)
            interval = RegionInterval(current_BB, start_stmt_idx, stmt_idx - 1)
            push!(intervals, interval)
            bb_interval[current_BB.cfgbb.label] = interval
            return true, interval
        end
    end
    interval = RegionInterval(current_BB, start_stmt_idx, last_stmt_idx) 
    push!(intervals, interval)
    bb_interval[current_BB.cfgbb.label] = interval
    
    # Current BB's statements have been scanned. Now successors
    for succ_BBcfg in current_BB.cfgbb.succs
        succ_BB = lives.basic_blocks[succ_BBcfg]
        if (DEBUG_LVL >= 2)
            println("looking at succ BB ", succ_BB.cfgbb.label, " whose dominators are ", loop_info.dom_dict[succ_BB.cfgbb.label])
        end
        
        if succ_BB.cfgbb.label == -2 || # the pseudo exit of the function
            (loop_bbs != nothing && !in(succ_BB.cfgbb.label, loop_bbs)) # out of loop TODO: allow region to include bbs out of loop. Just make sure it covers the whole loop, not part of it.
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, interval
            end
            exit = ExitEdge(current_BB, last_stmt_idx + 1, succ_BB, 0)
            push!(exits, exit)
            continue
        end            
        
        if !in(first_BB.cfgbb.label, loop_info.dom_dict[succ_BB.cfgbb.label])
            # first_BB does not dominate succ BB
            if loop_expected_on_every_path && !has_loop[current_BB.cfgbb.label]
                return false, interval
            end
            if in_loop[succ_BB.cfgbb.label]
                # we are going to insert reverse reodering on a new block between current and succ BB.
                # If succ BB is in a loop, the new block is also in a loop. We cannot affort that cost 
                return false, interval
            end
            
            exit = ExitEdge(current_BB, last_stmt_idx + 1, succ_BB, 0)
            push!(exits, exit)
            continue
        end
        
        if (visited[succ_BB.cfgbb.label])
            has_loop[current_BB.cfgbb.label] = has_loop[current_BB.cfgbb.label] || has_loop[succ_BB.cfgbb.label]
            push!(interval.succs, bb_interval[succ_BB.cfgbb.label])
            push!(bb_interval[succ_BB.cfgbb.label].preds, interval)
            continue
        end
        
        success, successor_interval = DFSGrowRegion(first_BB, succ_BB, 1, lives, in_loop, visited, has_loop, bb_interval, exits, intervals, symbolInfo, loop_info, loop_bbs)
        if !success
            return false, successor_interval
        else
            push!(interval.succs, successor_interval)
            push!(successor_interval.preds, interval)
        end
    end
    return true, interval
end

function growRegion(first_BB :: CompilerTools.CFGs.BasicBlock, 
                    mmread_stmt_idx :: Int64, 
                    lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                    in_loop :: Dict{Int,Bool}, 
                    symbolInfo :: Dict{Union(Symbol,Integer),Any}, 
                    loop_info :: CompilerTools.Loops.DomLoops, 
                    loop_bbs)
    visited  = Dict{Int, Bool}() # BB has been visited?
    has_loop = Dict{Int, Bool}() # Every path starting from the BB crosses a loop?
    bb_interval = Dict{Int, RegionInterval}()
    for (j,bb) in lives.cfg.basic_blocks 
        visited[bb.label]  = false
        has_loop[bb.label] = false
    end
    
    # Imagine there are always two invisible statements in any BB, whose statement indices are 0 and last_stmt_idx+1
    # So any real statement always has a statement before and after it. That is why we can do mmread_stmt_idx + 1 here
    # Similary, we can do any stmt_indx -1 as well.
    region = Region(lives.basic_blocks[first_BB], mmread_stmt_idx, Set{ExitEdge}(), RegionInterval[])
    success, interval = DFSGrowRegion(lives.basic_blocks[first_BB], lives.basic_blocks[first_BB], mmread_stmt_idx + 1, lives, in_loop, visited, has_loop, bb_interval, region.exits, region.intervals, symbolInfo, loop_info, loop_bbs)

    if (DEBUG_LVL >= 2)
        println("DFSGrowRegion successful?: ", success)
        println("Intervals of region found:")
        for interval in region.intervals
            show_interval(interval, 1)
            println("\tPreds: ")
            for pred in interval.preds
                show_interval(pred, 3)
            end
            println("\tSuccs: ")
            for succ in interval.succs
                show_interval(succ, 3)
            end
        end
    end
    
    if success
        return region, bb_interval
    else
        return nothing, bb_interval
    end
end

function regionFormationBasedOnMmread(lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                                      loop_info :: CompilerTools.Loops.DomLoops, 
                                      symbolInfo :: Dict{Union(Symbol,Integer),Any})
    in_loop = BBsInLoop(lives, loop_info)
    mmread_BB, mmread_stmt_idx = findMmread(lives, in_loop)
    if mmread_BB == nothing
        return nothing, nothing
    end
    return growRegion(mmread_BB, mmread_stmt_idx, lives, in_loop, symbolInfo, loop_info, nothing)
end

function regionFormationBasedOnLoop(L, lives, loop_info, symbolInfo)
    in_loop = BBsInLoop(lives, loop_info)
    return growRegion(lives.cfg.basic_blocks[L.head], 0, lives, in_loop, symbolInfo, loop_info, L.members)
end

# Old version, not differentiating LHS and RHS symbols. To remove
function findIAs(region::Region, symbolInfo)
    # Build inter-dependent arrays
    IAs = Dict{Any, Set}()
    for interval in region.intervals
        BB = interval.BB
        for stmt_index = interval.from_stmt_idx : interval.to_stmt_idx
            stmt = BB.statements[stmt_index]
            IAs[stmt] = Set{Any}()
            seeds = Set{Any}()
            push!(seeds, stmt.expr)
            while !isempty(seeds)
                seed = pop!(seeds)
                IA = Set{Any}()
                findIA(seed, IA, seeds, symbolInfo, 1)
                if !isempty(IA)
                    push!(IAs[stmt], IA)
                end
            end
#            println("###Stmt: ", stmt.expr)
#            println("  IAs:", IAs[stmt])        
        end
    end
    return IAs
end

# a set of inter-dependent arrays in the same statement. Some of them appear in LHS, some in RHS 
type InterDependentArrays
    LHS:: Set{Any}
    RHS:: Set{Any}
    InterDependentArrays() = new(Set(), Set())
end

function addToInterDependentArrays(currentNode, LHS::Bool, IA::InterDependentArrays, symbolInfo)
    if typeOfNode(currentNode, symbolInfo) <: AbstractArray
        if typeof(currentNode) == Symbol
            push!(LHS ? IA.LHS : IA.RHS, currentNode)
        elseif typeof(currentNode) == SymbolNode
            push!(LHS ? IA.LHS : IA.RHS, currentNode.name)
        end
    end 
end

function findInterDependentArrays(currentNode, LHS::Bool, IA::InterDependentArrays, seeds, symbolInfo, level)
#   printSpaces(level)
#   println("findIA: ", currentNode)
    addToInterDependentArrays(currentNode, LHS, IA, symbolInfo)
    if typeof(currentNode) <: Expr
        head = currentNode.head
        if head == :call || head == :call1
            head = currentNode.head
            args = currentNode.args

            # The first argument is the function, the others are the arguments for it.
            arg_types = ntuple(length(args) - 1, i-> typeOfNode(args[i+1], symbolInfo))
            all_numbers, some_arrays = analyze_types(currentNode.typ, arg_types)
            if all_numbers || !some_arrays
                # Result and args are all numbers, or there may be other types (like 
                # Range{UInt64}) but no regular arrays. No inter-dependent arrays.
                # Do nothing
            else
                module_name, function_name = resolve_module_function_names(args)
                if function_name == ""
                    if KEEP_GOING
                        show(UnresolvedFunction(head, args[1]))
                        @goto look_at_args
                    else
                        throw(UnresolvedFunction(head, args[1]))
                    end
                end
                fd = lookup_function_description(module_name, function_name, arg_types)
                if fd != nothing
                    for S in fd.IA
                        for x in S
                            if x == 0 # result of the function call
                                #findInterDependentArrays(currentNode, LHS, IA, seeds, symbolInfo, level+1)
                            else 
                                findInterDependentArrays(currentNode.args[x + 1], LHS, IA, seeds, symbolInfo, level+1)
                            end
                        end
                    end
                else
                    if KEEP_GOING
                        show(UndescribedFunction(module_name, function_name, arg_types))
                        @goto look_at_args
                    else
                        throw(UndescribedFunction(module_name, function_name, arg_types))
                    end
                end
            end
        end
        
@label look_at_args
        for c in currentNode.args
#            printSpaces(level+1)
#            println("c=", c, " type=", typeOfNode(c, symbolInfo))
            if typeOfNode(c, symbolInfo) <: Number
                push!(seeds, tuple(c, LHS))
            else
                findInterDependentArrays(c, LHS, IA, seeds, symbolInfo, level+1)
            end
        end
    end
end

function optimizeRegionCalls(region :: Region)
    for interval in region.intervals
        BB = interval.BB
        for stmt_index = 1 : length(BB.cfgbb.statements)
            BB.cfgbb.statements[stmt_index].expr = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(BB.cfgbb.statements[stmt_index].expr, optimize_calls, nothing))
        end
    end 
end

function findInterDependentArrays(region :: Region, symbolInfo :: Dict{Union(Symbol,Integer),Any})
    # Build inter-dependent arrays
    IAs = Dict{Any, Set}()
    for interval in region.intervals
        BB = interval.BB
        for stmt_index = interval.from_stmt_idx : interval.to_stmt_idx
            stmt = BB.statements[stmt_index]
            IAs[stmt] = Set{Any}()
            seeds = Set{Tuple{Any, Bool}}()
            push!(seeds, tuple(stmt.tls.expr, false))
            while !isempty(seeds)
                seed = pop!(seeds)
                IA = InterDependentArrays()
                if typeof(seed[1]) == Expr && seed[1].head == :(=)
                    # TODO: check that an assignment can happen only at top level, not in LHS or RHS 
                    findInterDependentArrays(seed[1].args[1], true, IA, seeds, symbolInfo, 1)
                    findInterDependentArrays(seed[1].args[2], false, IA, seeds, symbolInfo, 1)
                else
                    findInterDependentArrays(seed[1], seed[2], IA, seeds, symbolInfo, 1)
                end
                if !isempty(IA.LHS) || !isempty(IA.RHS)
                    push!(IAs[stmt], IA)
                end
            end
#            println("###Stmt: ", stmt.tls.expr)
#            println("  IAs:", IAs[stmt])        
        end
    end
    return IAs
end

function show_IAs(region::Region, IAs, prefix)
    println("******************************* ", prefix, " ***************************")
    for interval in region.intervals
        BB = interval.BB
        for stmt_index = interval.from_stmt_idx : interval.to_stmt_idx
            stmt = BB.statements[stmt_index]
            print(stmt)
            for IA in IAs[stmt]
                println("\tIA: L(", IA.LHS, ") R(", IA.RHS, ")")
            end
        end
    end
end

# A node in a graph for reorderable matrix discovery analysis (RMD)
# There are 4 kinds of nodes: 
# (1) Entry. It is inserted into the region as a special INSIDE statement.
# (2) Empty. It represents an empty interval INSIDE the region.
# (3) Normal. It represents a statement INSIDE the region.
# (4) Outside. It represents an OUTSIDE statement.
# The difference between INSIDE and OUTSIDE statement is that the INSIDE statement 
# is subject to the dataflow propagation, while the OUTSIDE statement is not.
# The OUTSIDE statement only provides a fixed IN (empty set) to its INSIDE predecessors
# and a fixed OUT (empty set) to its INSIDE successors.
const RMD_NODE_ENTRY = 0xfff0
const RMD_NODE_EMPTY = 0xfff1
const RMD_NODE_NORMAL = 0xfff2
const RMD_NODE_OUTSIDE = 0xfff4

FIRST_RMD_NODE_INDEX = 0
type RMDNode
    node_idx
    bbnum   
    stmt_idx
    kind      
    succs     :: Set{RMDNode}
    preds     :: Set{RMDNode}
    In        :: Set{Any} #Set{Symbol} causes incompatibility with TopLevelStatement.def and use TODO: ask Todd to fix TopLevelStatement declaration
    Out       :: Set{Any} #TODO: change back to Set{Symbol} 
            
    RMDNode(basicBlock_num, statement_idx, node_kind) = (
        global FIRST_RMD_NODE_INDEX;
        FIRST_RMD_NODE_INDEX += 1;
        new(FIRST_RMD_NODE_INDEX, basicBlock_num, statement_idx, node_kind,
            Set{RMDNode}(), Set{RMDNode}(), Set{Any}(), Set{Any}()))
end

function show_RMD_node_kind(node::RMDNode)
        if node.kind == RMD_NODE_ENTRY
            print(" ENTRY ")
        elseif node.kind == RMD_NODE_EMPTY
            print(" EMPTY ")
        elseif node.kind == RMD_NODE_NORMAL
            print(" NORMAL ")
        elseif node.kind == RMD_NODE_OUTSIDE
            print(" OUTSIDE ")            
        else
            print(" WRONG KIND!!!! ")
        end
end

function show_RMD_node(node::RMDNode, lives)
        print(node.node_idx, " <BB ", node.bbnum, ", Stmt ", node.stmt_idx, ">")
        show_RMD_node_kind(node)
        print(" Preds(")
        for pred in node.preds
            print(pred.node_idx, " ")
        end
        print(")")
            
        print(" Succs(")
        for succ in node.succs
            print(succ.node_idx, " ")
        end
        print(")")
        
        print("  In(")
        for item in node.In
            print(item, " ")
        end
        print(")")

        print("  Out(")
        for item in node.Out
            print(item, " ")
        end
        println(")")
        
        if node.kind == RMD_NODE_NORMAL
            BB = lives.basic_blocks[lives.cfg.basic_blocks[node.bbnum]]
            stmt = BB.statements[node.stmt_idx]
            println("\t", stmt, " ", stmt.tls)
        end
end

function show_RMD_graph(nodes, outside_nodes, lives, prefix::String)
    println("******************************* ", prefix, " ***************************")
    for node in union(nodes, outside_nodes)
        show_RMD_node(node, lives)
    end
end

const UNIVERSE_SYM = gensym("universe")

function intersect_preds_out(node)
    S = Set{Any}()
    first = true
    for pred in node.preds
        if first || in(UNIVERSE_SYM, S)
            S = copy(pred.Out)
            first = false
        else
            if !in(UNIVERSE_SYM, pred.Out)
                S = intersect(S, pred.Out)
            end 
        end
    end
    S
end

function intersect_succs_in(node)
    S = Set{Any}()
    first = true
    for succ in node.succs
        if first || in(UNIVERSE_SYM, S)
            S = copy(succ.In)
            first = false
        else
            if !in(UNIVERSE_SYM, succ.In)
                S = intersect(S, succ.In)
            end 
        end
    end
    S
end

function forward_transfer(B::RMDNode, S::Set, IAs, lives)
    assert(B.kind != RMD_NODE_OUTSIDE)
    if B.kind == RMD_NODE_ENTRY || B.kind == RMD_NODE_EMPTY
        return copy(S)
    end
    
    BB = lives.basic_blocks[lives.cfg.basic_blocks[B.bbnum]]
    stmt = BB.statements[B.stmt_idx]
    if isempty(IAs[stmt])
        return copy(S)
    end
    
    result = Set()
    for x in S
        for IA in IAs[stmt]
            if in(x, IA.RHS)
                union!(result, IA.RHS)
                # no matter a RHS symbol is in LHS or not, it should be transfered:
                # If it is in, then of course it will (It means the new def of the symbol)
                # If not, it will as well (It means the current def of the symbol)
                union!(result, IA.LHS)
            else
                if in(x, IA.LHS)
                    # The current def of the symbol is killed. The new def, unless
                    # some symbol in the RHS is inter-dependent on it, will not be added
                    # Do nothing.
                else
                    push!(result, x)
                end
            end
        end
    end
    return result
end

function backward_transfer(B::RMDNode, S::Set, IAs, lives)
    assert(B.kind != RMD_NODE_OUTSIDE)
    if B.kind == RMD_NODE_ENTRY || B.kind == RMD_NODE_EMPTY
        return copy(S)
    end
    
    BB = lives.basic_blocks[lives.cfg.basic_blocks[B.bbnum]]
    stmt = BB.statements[B.stmt_idx]
    if isempty(IAs[stmt])
        return copy(S)
    end
    
    result = Set()
    for x in S
        for IA in IAs[stmt]
            if in(x, IA.LHS)
                union!(result, IA.RHS)
            else
                if in(x, IA.RHS)
                    union!(result, IA.RHS)
                else
                    push!(result, x)
                end
            end
        end
    end
    return result
end

function forwardPass(nodes, IAs, initialization::Bool, lives)
    ever_changed = false
    changed = true
    while changed
        changed = false
        for node in nodes
            if initialization && node.kind == RMD_NODE_ENTRY
                continue
            end
            S = intersect_preds_out(node)
            if !initialization
                S = union(node.In, S)
            end
            if !(S == node.In)
                changed = true
                ever_changed = true
                node.In = S
            end
            S = forward_transfer(node, node.In, IAs, lives)
            if !initialization
                S = union(node.Out, S)
            end
            if !(S == node.Out)
                changed = true
                ever_changed = true
                node.Out = S
            end
        end
    end
    ever_changed
end

function backwardPass(nodes, IAs, lives)
    ever_changed = false
    changed = true
    while changed
        changed = false
        for node in nodes
            S = intersect_succs_in(node)
            S = union(node.Out, S)
            if !(S == node.Out)
                changed = true
                ever_changed = true
                node.Out = S
            end
            S = backward_transfer(node, node.Out, IAs, lives)
            S = union(node.In, S)
            if !(S == node.In)
                changed = true
                ever_changed = true
                node.In = S
            end
        end
    end
    ever_changed
end

function intervalsAreAdjacent(pred_interval::RegionInterval, interval::RegionInterval)
    assert(interval.from_stmt_idx >= 1)
    if interval.from_stmt_idx == 1 
        assert(pred_interval.to_stmt_idx <= length(pred_interval.BB.statements))
        if pred_interval.to_stmt_idx == length(pred_interval.BB.statements)
            return true # true even if pred_interval is empty (from_stmt_idx > to_stmt_idx)
        else 
            return false
        end
    else
        return false
    end
end

function makeConnectPredOutsideNode(interval :: RegionInterval, predBB :: CompilerTools.LivenessAnalysis.BasicBlock, pred_stmt_idx, first_node, outside_nodes)
    outside_node = RMDNode(predBB.cfgbb.label, pred_stmt_idx, RMD_NODE_OUTSIDE)
    push!(outside_nodes, outside_node)
    push!(first_node[interval].preds, outside_node)
    push!(outside_node.succs, first_node[interval])
end

function makeConnectSuccOutsideNode(interval::RegionInterval, succBB :: CompilerTools.LivenessAnalysis.BasicBlock, succ_stmt_idx, last_node, outside_nodes)
    outside_node = RMDNode(succBB.cfgbb.label, succ_stmt_idx, RMD_NODE_OUTSIDE)
    push!(outside_nodes, outside_node)
    push!(last_node[interval].succs, outside_node)
    push!(outside_node.preds, last_node[interval])
end

function reorderableMatrixDiscovery(
                       lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                       loop_info :: CompilerTools.Loops.DomLoops, 
                       region :: Region, 
                       bb_interval :: Dict{Int, RegionInterval}, 
                       IAs :: Dict{Any,Set},
                       root)
    if isempty(region.intervals) 
        return nothing, nothing, nothing
    end
    
    # build a graph, where each node is a statement.
    first_node = Dict{RegionInterval, RMDNode}()
    last_node  = Dict{RegionInterval, RMDNode}()
    nodes      = RMDNode[] # All INSIDE nodes

    # build a pseudo entry and exit node
    entry = RMDNode(0, 0, RMD_NODE_ENTRY)
    push!(nodes, entry)

    # TODO visit the intervals in reverse postorder, so that the nodes
    # can be processed in a natural order
    for interval in region.intervals
        if interval.from_stmt_idx > interval.to_stmt_idx
            # Empty interval. Make a pseudo node
            node = RMDNode(interval.BB.cfgbb.label, interval.from_stmt_idx, RMD_NODE_EMPTY)
            push!(nodes, node)
            first_node[interval] = last_node[interval] = node
        else
            prev = nothing    
            for stmt_idx in interval.from_stmt_idx : interval.to_stmt_idx
                node = RMDNode(interval.BB.cfgbb.label, stmt_idx, RMD_NODE_NORMAL)
                push!(nodes, node)
                
                if stmt_idx == interval.from_stmt_idx
                    first_node[interval] = node
                end
                if stmt_idx == interval.to_stmt_idx
                    last_node[interval] = node
                end
                
                if prev != nothing
                    push!(prev.succs, node)
                    push!(node.preds, prev)
                end
                
                prev = node
            end
        end
    end

    # Now all INSIDE nodes have been created, and those in the same interval have been
    # connected. Connected nodes between intervals, and create OUTSIDE predecessor and
    # successors
    outside_nodes = RMDNode[]
    first_interval = region.intervals[1]
    push!(entry.succs, first_node[first_interval])
    push!(first_node[first_interval].preds, entry)

    for interval in region.intervals
        BB = interval.BB
        if interval.from_stmt_idx == 1
            # interval covers the start of the BB. Connect to predecessors
            # Exception: if this is a loop region and BB is the loop head, then ENTRY has
            # represented all its predecessors except the backedge one. So we need only to
            # build the predecessor on the backedge.
            BB_is_loop_head = ((region.first_BB == BB) && (region.mmread_stmt_idx == 0))
            for predcfg in BB.cfgbb.preds
                pred = lives.basic_blocks[predcfg]
                if BB_is_loop_head && !in(BB.cfgbb.label, loop_info.dom_dict[pred.cfgbb.label])
                    continue
                end
                if haskey(bb_interval, pred.cfgbb.label)
                    pred_interval = bb_interval[pred.cfgbb.label]
                    if intervalsAreAdjacent(pred_interval, interval)
                        push!(last_node[pred_interval].succs, first_node[interval])
                        push!(first_node[interval].preds, last_node[pred_interval])
                    else 
                        makeConnectPredOutsideNode(interval, pred, length(pred.statements), first_node, outside_nodes)
                    end
                else
                    # pred has no interval. Make an outside node for it.
                    makeConnectPredOutsideNode(interval, pred, length(pred.statements), first_node, outside_nodes)
                end
            end
        else
            assert(interval.from_stmt_idx > 1)
            # interval starts from the middle of the BB. Make and connect to an OUTSIDE predecessor in the current BB
            makeConnectPredOutsideNode(interval, BB, interval.from_stmt_idx - 1, first_node, outside_nodes)
        end
        
        if interval.to_stmt_idx == length(BB.statements)
            # interval covers the end of the BB. Connect to successors
            for succcfg in BB.cfgbb.succs
                succ = lives.basic_blocks[succcfg]
                if haskey(bb_interval, succ.cfgbb.label)
                    succ_interval = bb_interval[succ.cfgbb.label]
                    if intervalsAreAdjacent(interval, succ_interval)
                        push!(last_node[interval].succs, first_node[succ_interval])
                        push!(first_node[succ_interval].preds, last_node[interval])
                    else 
                        makeConnectSuccOutsideNode(interval, succ, 1, last_node, outside_nodes)
                    end
                else
                    # succ has no interval. Make an outside node for it.
                    makeConnectSuccOutsideNode(interval, succ, 1, last_node, outside_nodes)
                end
            end
        else
            assert(interval.to_stmt_idx < length(BB.statements))
            # interval stops at the middle of the BB. Make and connect to an OUTSIDE successor in the current BB
            makeConnectSuccOutsideNode(interval, BB, interval.to_stmt_idx + 1, last_node, outside_nodes)
        end
    end
    
    # The region must have at least one exit to outside.
    assert(!isempty(outside_nodes))
    
    # Do bidirectional dataflow analysis on the graph
    if (DEBUG_LVL >= 2)
        show_RMD_graph(nodes, outside_nodes, lives, "Initial RMD graph:")
    end

    # Step 1: initialization
    for node in nodes
        if node.kind == RMD_NODE_ENTRY
            entry.In  = Set{Any}()
            entry.Out = Set{Any}()
            push!(entry.In, root)
            push!(entry.Out, root)
         else
            node.Out = Set{Any}()
            push!(node.Out, UNIVERSE_SYM)
        end
    end

    for outside_node in outside_nodes
        outside_node.In  = Set{Any}()
        outside_node.Out = Set{Any}()
    end

    forwardPass(nodes, IAs, true, lives)
        if (DEBUG_LVL >= 2)
            show_RMD_graph(nodes, outside_nodes, lives, "RMD graph after initialization:")
        end
    
    # repetitive backward and forward pass
    changed = true
    while changed
        changed = backwardPass(nodes, IAs, lives)
        changed |= forwardPass(nodes, IAs, false, lives)
        if (DEBUG_LVL >= 2)
            show_RMD_graph(nodes, outside_nodes, lives, "RMD graph after 1 iteration:")
        end
    end
    assert(entry.In == entry.Out)
        
    return nodes, entry, outside_nodes
end

function insertStatementsOnEdge(lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                                new_stmts, 
                                node::RMDNode, 
                                node_BB :: CompilerTools.LivenessAnalysis.BasicBlock, 
                                succ::RMDNode, 
                                succ_BB :: CompilerTools.LivenessAnalysis.BasicBlock)
    assert(node.kind != RMD_NODE_ENTRY) # Reordering from entry to the first BB has been processed specially. Not here.

    if isempty(new_stmts)
        return
    end

    BB = nothing
    insert_at = 0
    new_goto_stmt = nothing
    if node_BB == succ_BB
        #insert inside the same block
        assert(node.stmt_idx + 1 == succ.stmt_idx)
        assert(succ.stmt_idx >= 1)
        assert(succ.stmt_idx <= length(succ_BB.statements) + 1)
        BB = node_BB.cfgbb
        insert_at = succ.stmt_idx
    else
        (BB, new_goto_stmt) = CompilerTools.CFGs.insertBetween(lives.cfg, node.bbnum, succ.bbnum)
        lives.basic_blocks[BB] = CompilerTools.LivenessAnalysis.BasicBlock(BB)
        insert_at = 1
    end
    
    assert(typeof(BB) == CompilerTools.CFGs.BasicBlock)

    i = 0
    for new_stmt in new_stmts
        insert!(BB.statements, insert_at + i, CompilerTools.CFGs.TopLevelStatement(0, new_stmt))
        #insert!(BB.statements, insert_at + i, CompilerTools.LivenessAnalysis.TopLevelStatement(0, new_stmt))
        #CompilerTools.CFGs.insertat!(BB.statements, CompilerTools.CFGs.TopLevelStatement(0, new_stmt), insert_at + i)
        i += 1
    end
    if new_goto_stmt != nothing
        push!(BB.statements, new_goto_stmt)
    end
end

function nodeLiveIn(node::RMDNode, BB :: CompilerTools.LivenessAnalysis.BasicBlock)
    if node.stmt_idx < 1
        assert(node.stmt_idx == 0)
        return BB.live_in
    elseif node.stmt_idx > length(BB.statements)
        assert(node.stmt_idx == (length(BB.statements) + 1))
        return BB.live_out
    else
        stmt = BB.statements[node.stmt_idx]
        return stmt.live_in
    end
end

function regionTransformation(funcAST, 
                       lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                       loop_info :: CompilerTools.Loops.DomLoops, 
                       symbolInfo :: Dict{Union(Symbol,Integer),Any}, 
                       region :: Region, 
                       bb_interval :: Dict{Int, RegionInterval}, 
                       IAs :: Dict{Any,Set},
                       M, 
                       back_edge)
    # The region is based on either a mmread(), in which case LHS of the statement will be the root symbol, 
    # or a loop, in which case the root symbol is identified from the function's argument list
    assert(region.mmread_stmt_idx != 0 || M != nothing)
    
    if region.mmread_stmt_idx != 0
        mmread_stmt = region.first_BB.statements[ region.mmread_stmt_idx ] 
        M = mmread_stmt.tls.expr.args[1]
    end 
        
    nodes, entry, outside_nodes = reorderableMatrixDiscovery(lives, loop_info, region, bb_interval, IAs, M)
    
    #TODO: benefit-cost analysis
    
    # New a vector to hold the new reordering statements R(LiveIn) before the region. 
    new_stmts_before_region = Expr[]

    if region.mmread_stmt_idx != 0
        reorderedBeforeRegion = intersect(entry.Out, mmread_stmt.live_out)
        dprintln(2, "To be reordered before region: ", reorderedBeforeRegion)

        if isempty(reorderedBeforeRegion)
            return
        end
        
        # Replace the original "M=mmread(filename)" with
        #   T   = mmread_reorder(filename) // T is a tuple
        #   M = (top(tupleref))(T, 1)
        #   P   = (top(tupleref))(T, 2)
        #   P'  = (top(tupleref))(T, 3)
        T = gensym("T")
        P = gensym("P")
        Pprime = gensym("Pprime")
        
        # TODO: change the type of the original RHS to the right type
        mmread_stmt.tls.expr.args[1] = T
        mmread_stmt.tls.expr.args[2].args[1].args[3] = QuoteNode(:mmread_reorder) # change mmread to mmread_reorder
        
        stmt = Expr(:(=), M, Expr(:call, TopNode(:tupleref), T, 1))
        push!(new_stmts_before_region, stmt)
        
        stmt = Expr(:(=), P, Expr(:call, TopNode(:tupleref), T, 2))
        push!(new_stmts_before_region, stmt)
        
        stmt = Expr(:(=), Pprime, Expr(:call, TopNode(:tupleref), T, 3))
        push!(new_stmts_before_region, stmt)
    else
        # The symbols that should be reordered before the region are the reordered uses live into the region
        reorderedBeforeRegion = intersect(entry.Out, region.first_BB.live_in)
        dprintln(2, "To be reordered before region: ", reorderedBeforeRegion)

        if isempty(reorderedBeforeRegion)
            return
        end
    
        # Allocate space to store the permutation and inverse permutation info
        (P, Pprime) = allocateForPermutation(M, new_stmts_before_region)
    
        # Compute P and Pprime, and reorder M
        reorderMatrix(M, P, Pprime, new_stmts_before_region, true, true, true)
    
    end
    
    # Now reorder other arrays
    for sym in reorderedBeforeRegion
        if sym != M
            if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                reorderMatrix(sym, P, Pprime, new_stmts_before_region, false, true, true)
            else
                reorderVector(sym, P, new_stmts_before_region)
            end
        end
    end

    if region.mmread_stmt_idx != 0
        i = 1
        for new_stmt in new_stmts_before_region
            insert!(region.first_BB.cfgbb.statements, region.mmread_stmt_idx + i, CompilerTools.CFGs.TopLevelStatement(0, new_stmt))    
            i += 1
        end
    else 
        (new_bb, new_goto_stmt) = CompilerTools.CFGs.insertBefore(lives.cfg, region.first_BB.cfgbb.label, true, back_edge)
        lives.basic_blocks[new_bb] = CompilerTools.LivenessAnalysis.BasicBlock(new_bb)
        for new_stmt in new_stmts_before_region
            CompilerTools.CFGs.addStatementToEndOfBlock(lives.cfg, new_bb, new_stmt)
        end
        if new_goto_stmt != nothing
          push!(new_bb.statements, new_goto_stmt)
        end
    end

    # Only those arrays that are updated in the region need to be reverse reordered afterwards
    updatedInRegion = Set{Symbol}()
    for interval in region.intervals
        BB = interval.BB
        for stmt_index = interval.from_stmt_idx : interval.to_stmt_idx
            stmt = BB.statements[stmt_index]
            union!(updatedInRegion, stmt.def)
        end
    end
    dprintln(2, "updatedInRegion: ", updatedInRegion)
 
    # Now process other nodes. Consider outside nodes as well: we may need to insert on
    # an edge between an inside and outside node.
    all_nodes = union(nodes, outside_nodes)
    # Look at all edges in the graph
    for node in all_nodes
        if node == entry # Entry has been processed before
            continue
        end
        # an inside node must have at least one successor
        assert(node.kind == RMD_NODE_OUTSIDE || !isempty(node.succs))
        BB = lives.basic_blocks[lives.cfg.basic_blocks[node.bbnum]]
        for succ in node.succs
            succ_BB = lives.basic_blocks[lives.cfg.basic_blocks[succ.bbnum]]
            succ_live_in = nodeLiveIn(succ, succ_BB)
            
            # compute what to be reorder on this edge
            reorder = setdiff(intersect(succ.In, succ_live_in), node.Out)
            new_stmts = Expr[]
            for sym in reorder
                if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                    reorderMatrix(sym, P, Pprime, new_stmts, false, true, true)
                else
                    reorderVector(sym, P, new_stmts)
                end
            end
            insertStatementsOnEdge(lives, new_stmts, node, BB, succ, succ_BB)

            # compute what to be reverse reordered on this edge
            reverseReorder = setdiff(intersect(node.Out, succ_live_in, updatedInRegion), succ.In)
            new_stmts = Expr[]
            for sym in reverseReorder
                if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                    reverseReorderMatrix(sym, P, Pprime, new_stmts, false, true, true)
                else
                    reverseReorderVector(sym, P, new_stmts)
                end
            end
            insertStatementsOnEdge(lives, new_stmts, node, BB, succ, succ_BB)
        end
    end

    if (DEBUG_LVL >= 2)
        println("******** CFG after region transformation: ********")
        show(lives.cfg);
    end 
end

function reorderRegion(funcAST, 
                       lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                       loop_info :: CompilerTools.Loops.DomLoops, 
                       symbolInfo :: Dict{Union(Symbol,Integer),Any}, 
                       region :: Region, 
                       bb_interval :: Dict{Int, RegionInterval}, 
                       M, 
                       back_edge)
    IAs = findInterDependentArrays(region, symbolInfo)
    if (DEBUG_LVL >= 2)
        show_IAs(region, IAs, "Inter-dependent arrays")
    end
    regionTransformation(funcAST, lives, loop_info, symbolInfo, region, bb_interval, IAs, M, back_edge)
    optimizeRegionCalls(region)
end

function reorder(funcAST, 
                 lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                 loop_info :: CompilerTools.Loops.DomLoops, 
                 symbolInfo :: Dict{Union(Symbol,Integer),Any})
    if (DEBUG_LVL >= 2)
        println("******** CFG before reorder: ********")
        show(lives.cfg);
    end 
    
    assert(funcAST.head == :lambda)
    args = funcAST.args
    assert(length(args) == 3)
    assert(typeof(args[3]) == Expr && args[3].head == :body)

    # First, see if there is an mmread() outside all loops.
    # If not, select a sparse matrix from the function AST's arguments. 
    # So far, choose the first sparse matrix in the arguments.
    # TODO: have a heuristic to choose the best candidate? Or explicitly
    # point out the candidate with user annotation?
    # TODO: it really does not matter whether the sparse matrix is from
    # an mmread or from an argument -- from the perspective of region formation,
    # identifying IAs, and region transformation. So should unify the two cases.
    region, bb_interval = regionFormationBasedOnMmread(lives, loop_info, symbolInfo)
    if region != nothing 
        reorderRegion(funcAST, lives, loop_info, symbolInfo, region, bb_interval, nothing, nothing)    
        return
    end

    local param = args[1]
    found = false
    M = nothing
    for i = 1:length(param)
        if typeOfNode(param[i], symbolInfo) <: AbstractSparseMatrix
            found = true
            M = param[i]
            break
        end
    end    
    if found 
        # form a region for each outermost loop in the function
        in_loop = BBsInLoop(lives, loop_info)
        for L in loop_info.loops
            is_outermost = true
            for L1 in loop_info.loops
                assert(L == L1 || L.members != L1.members)
                if L != L && L.members < L1.members
                    is_outermost = false
                    break
                end
            end
            if is_outermost
                region, bb_interval = regionFormationBasedOnLoop(L, lives, loop_info, symbolInfo)
                if region != nothing
                    reorderRegion(funcAST, lives, loop_info, symbolInfo, region, bb_interval, M, L.back_edge)
                end
            end
        end
    end
end

type CallSite
    args # args of the call
    matrices :: Set # matrices in the arguments of the call. They are symbols or GenSym's
    fknob_creator # A library call to create a function knob for this call site
    fknob_deletor # A library call to delete the function knob for this call site
end

type CallSites
    sites :: Set{CallSite}
    symbolInfo :: Dict
end

function context_may_help(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    if typeof(ast) <: Expr
        head = ast.head
        if head == :call || head == :call1
            args = ast.args
            module_name, function_name = resolve_module_function_names(args)
            if module_name == "" && function_name == "fwdTriSolve!" &&
                length(args) == 3 && 
                typeOfNode(args[2], call_sites.symbolInfo) == SparseMatrixCSC &&
                typeOfNode(args[3], call_sites.symbolInfo) == Vector
                    fknob_creator = (:NewForwardTriangularSolveKnob, LIB_PATH)
                    fknob_deletor = (:DeleteForwardTriangularSolveKnob, LIB_PATH)
                    site::CallSite = CallSite(args, Set(), fknob_creator, fknob_deletor)
                    if typeof(args[2]) == SymbolNode
                        push!(site.matrices, args[2].name)
                    else
                        push!(site.matrices, args[2])
                    end
                    push!(call_sites, site) 
            end
        end
    end
    return nothing
end

function index_of_stmt(BB, stmt)
    for index = 1 : length(BB.cfgbb.statements)
        if BB.cfgbb.statements[index] == stmt
            return index
        end
    end
    return -1
end

function insert_new_stmt_before(new_stmt, BB, stmt)
    index = index_of_stmt(BB, stmt)
    assert(index > 0)
    insert!(BB.cfgbb.statements, index, new_stmt)
end

function insert_new_matrix_knob(M, func_entry_BB, insert_at)
    mknob = gensym(string("mknob", string(M)))
    new_stmt = Expr(:(=), mknob,
                Expr(:call, GlobalRef(SparseAccelerator, :new_matrix_knob)))
    insert!(func_entry_BB.statements, insert_at, CompilerTools.CFGs.TopLevelStatement(0, new_stmt))
    mknob
end

function insert_increment_matrix_version(mknob, stmt, stmt_BB)
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :increment_matrix_version), mknob)
    insert_new_stmt_before(new_stmt, stmt_BB, stmt)
end

function insert_new_function_knob(call_site, func_entry_BB, matrix_knobs, insert_at)
    fknob = GenSym("fknob")
    new_stmt = Expr(:(=), fknob, 
                Expr(:call, GlobalRef(SparseAccelerator, :new_function_knob), 
                    call_site.fknob_creator)
               )
    insert!(func_entry_BB.statements, insert_at, CompilerTools.CFGs.TopLevelStatement(0, new_stmt))
    push!(call_site.args, fknob)
    fknob
end

function delete_function_knob(fknob_deletor, fknob, exit_BB, exit_stmt)
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :delete_function_knob), fknob_deletor, fknob)
    insert_new_stmt_before(new_stmt, exit_BB, exit_stmt)
end

function delete_matrix_knob(mknob, exit_BB, exit_stmt)
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :delete_matrix_knob), mknob)
    insert_new_stmt_before(new_stmt, exit_BB, exit_stmt)
end

# Analyze the sparse library function calls in the AST of a function, 
# and insert knobs(context-sensitive info)
function insert_knobs(funcAST, 
                      lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                      loop_info :: CompilerTools.Loops.DomLoops,
                      symbolInfo :: Dict{Union(Symbol,Integer),Any})
    if (DEBUG_LVL >= 2)
        println("******** CFG before insert_knobs: ********")
        show(lives.cfg);
    end 
    
    assert(funcAST.head == :lambda)
    args = funcAST.args
    assert(length(args) == 3)
    assert(typeof(args[3]) == Expr && args[3].head == :body)

    body = funcAST.args[3]

    # Scan all the matrices' defs and uses statements, and build up knobs for 
    # those related with routines that may take advantage of the context info,
    # like triangular solver, etc. 
    defs = Dict{Any, Set} # Map from a variable to a set of statements defining it
    uses = Dict{Any, Set{Tuple}} # Map from a variable to a set of (statement, args) tuples using it
    matrices = Set() # The matrices (variables) that should have context info
    call_sites = CallSites(Set{CallSite}(), symbolInfo)
    BB_in_loop, stmt_BB = statementsInLoop(lives, loop_info)
    exits = Set()
    for stmt in body.args
        for d in stmt.def
            if typeOfNode(d, symbolInfo) <: AbstractSparseMatrix
                push!(defs[d], stmt)
            end
        end
        
        # Only care about the reference points inside a loop: only there, library
        # can show benefit of reusing (every iteration)
        if BB_in_loop(stmt_BB[stmt]) # the stmt' BB is in a loop
            CompilerTools.AstWalker.AstWalk(stmt, context_may_help, call_sites)
        end
        
        expr = stmt.expr
        if typeof(expr) == Expr && (expr.head == :return || expr.head == :throw)
            push!(exits, stmt)
        end
    end

    func_entry_BB = lives.basic_blocks[-1]
    insert_at = 1
    matrix_knobs = Dict()
    for call_site in call_sites
        for M in call_site.matrices
            if !haskey(matrix_knobs, M)
                # Insert knob initialization at the beginning of the function
                mknob = insert_new_matrix_knob(M, func_entry_BB, insert_at)
                insert_at = insert_at + 1
                matrix_knobs[M] = mknob

                # Insert knob update before every statement that defines the matrix
                for stmt in defs[M]
                    insert_increment_matrix_version(mknob, stmt, stmt_BB)
                end
            end
        end
    end
    
    function_knobs = Set()
    for call_site in call_sites
        fknob = insert_new_function_knob(call_site, func_entry_BB, matrix_knobs, insert_at)
        insert_at = insert_at + 1
        push!(function_knobs, (fknob, call_site.fnob_deletor))
        for M in call_site.matrices
            add_matrix_knob(fknob, matrix_knobs[M])
        end
    end
    
    # Delete all the knobs at each exit of the function
    for exit_stmt in exits
        exit_BB = stmt_BB[exit_stmt]
        for (fknob, fnob_deletor) in function_knobs
            delete_function_knob(fnob_deletor, fknob, exit_BB, exit_stmt)
        end
        for mknob in values(matrix_knobs)
            delete_matrix_knob(mknob, exit_BB, exit_stmt)
        end
    end
    
    funcAST
end

# Try reordering, then insert knobs.
# TODO: combine the two together so that we scan AST only once.
function sparse_analyze(ast, 
                        lives :: CompilerTools.LivenessAnalysis.BlockLiveness, 
                        loop_info :: CompilerTools.Loops.DomLoops, 
                        symbolInfo :: Dict{Union(Symbol,Integer),Any})
    dprintln(2, "***************** Sparse analyze *****************")

    reorder(ast, lives, loop_info, symbolInfo)
    insert_knobs(ast, lives, loop_info, symbolInfo)

    body_reconstructed   = CompilerTools.CFGs.createFunctionBody(lives.cfg)
    ast.args[3].args = body_reconstructed

    flush(STDOUT::IO)
    
    return ast
end
