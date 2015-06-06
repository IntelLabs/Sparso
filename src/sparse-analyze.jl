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

const LIB_PATH = "../lib/libcsr.so"

# In reordering, we insert some calls to the following 3 functions. So they are executed secretly
# Reorder sparse matrix A and store the result in newA. A itself is not changed.
function CSR_ReorderMatrix(A::SparseMatrixCSC, newA::SparseMatrixCSC, P::Vector, Pprime::Vector, getPermutation::Bool, oneBasedInput::Bool, oneBasedOutput::Bool)
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
    elseif ast.args[1] == :(+)
      if length(ast.args) == 3
        dprintln(3,"optimize_calls found +")
        ast.args[2] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[2], optimize_calls, nothing))
        ast.args[3] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[3], optimize_calls, nothing))
        arg1 = ast.args[2]
        arg2 = ast.args[3]
        dprintln(3,"arg1 = ", arg1, " arg2 = ", arg2, " arg1.typ = ", deepType(arg1), " arg2.typ = ", deepType(arg2))
        if deepType(arg1) <: Vector && deepType(arg2) <: Vector 
          (is_mul1, scalar1, vec1) = isMul(arg1)
          (is_mul2, scalar2, vec2) = isMul(arg2)
          dprintln(3,"optimize_calls found vector/vector +. mul1 = ", is_mul1, " mul2 = ", is_mul2)
          if is_mul1 || is_mul2 || true
            orig_args = ast.args
            ast.args = Array(Any,5)
            dprintln(3,"optimize_calls converting to SparseAccelerator.WAXPBY")
            ast.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY))
            if is_mul1
              ast.args[2] = scalar1
              ast.args[3] = vec1
            else
              ast.args[2] = 1
              ast.args[3] = arg1
            end
            if is_mul2
              ast.args[4] = scalar2
              ast.args[5] = vec2
            else
              ast.args[4] = 1
              ast.args[5] = arg2
            end
          end
        end
        return ast
      end
    elseif ast.args[1] == :(-)
      if length(ast.args) == 3
        dprintln(3,"optimize_calls found -")
        ast.args[2] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[2], optimize_calls, nothing))
        ast.args[3] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[3], optimize_calls, nothing))
        arg1 = ast.args[2]
        arg2 = ast.args[3]
        dprintln(3,"arg1 = ", arg1, " arg2 = ", arg2, " arg1.typ = ", deepType(arg1), " arg2.typ = ", deepType(arg2))
        if deepType(arg1) <: Vector && deepType(arg2) <: Vector 
          (is_mul1, scalar1, vec1) = isMul(arg1)
          (is_mul2, scalar2, vec2) = isMul(arg2)
          dprintln(3,"optimize_calls found vector/vector +. mul1 = ", is_mul1, " mul2 = ", is_mul2)
          if is_mul1 || is_mul2 || true
            orig_args = ast.args
            ast.args = Array(Any,5)
            dprintln(3,"optimize_calls converting to SparseAccelerator.WAXPBY")
            ast.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY))
            if is_mul1
              ast.args[2] = scalar1
              ast.args[3] = vec1
            else
              ast.args[2] = 1
              ast.args[3] = arg1
            end
            if is_mul2
              ast.args[4] = CompilerTools.LivenessAnalysis.TypedExpr(deepType(scalar2), :call, :(-), scalar2)
              ast.args[5] = vec2
            else
              ast.args[4] = -1
              ast.args[5] = arg2
            end
          end
        end
        return ast
      end
    elseif ast.args[1] == :(*)
      dprintln(3,"optimize_calls found *")
      ast.args[2] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[2], optimize_calls, nothing))
      ast.args[3] = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(ast.args[3], optimize_calls, nothing))
      arg1 = ast.args[2]
      arg2 = ast.args[3]
      dprintln(3,"arg1 = ", arg1, " arg2 = ", arg2, " arg1.typ = ", deepType(arg1), " arg2.typ = ", deepType(arg2))
      if deepType(arg1) <: SparseMatrixCSC && deepType(arg2) <: Vector
        dprintln(3,"optimize_calls converting to SparseAccelerator.SpMV")
        ast.args[1] = CompilerTools.LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV))
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
        show(lives);
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
        for stmt_index = 1:length(bbs[bbnum].statements)
            # Replace calls to optimized versions provided in SparseAccelerator module.
            bbs[bbnum].statements[stmt_index].expr = CompilerTools.AstWalker.get_one(CompilerTools.AstWalker.AstWalk(bbs[bbnum].statements[stmt_index].expr, optimize_calls, nothing))
            stmt = bbs[bbnum].statements[stmt_index]
            IAs[stmt] = Set{Any}()
            seeds = Set{Any} ()
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
            for stmt in bbs[bbnum].statements
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
    (new_bb, new_goto_stmt) = CompilerTools.LivenessAnalysis.insertBefore(lives, L.head, true, L.back_edge)
    for stmt in new_stmts_before_L
        CompilerTools.LivenessAnalysis.addStatementToEndOfBlock(lives, new_bb, stmt)
    end
    if new_goto_stmt != nothing
      push!(new_bb.statements, new_goto_stmt)
    end
    
    for (pred, succ, new_stmts) in new_stmts_after_L
        (new_bb, new_goto_stmt) = CompilerTools.LivenessAnalysis.insertBetween(lives, pred, succ)
        for stmt in new_stmts
            CompilerTools.LivenessAnalysis.addStatementToEndOfBlock(lives, new_bb, stmt)
        end
        if new_goto_stmt != nothing
          push!(new_bb.statements, new_goto_stmt)
        end
    end

    if(DEBUG_LVL >= 2)
        println("******** CFG after reordering: ********")
        show(lives);
    end 
end

# Map each BB to a boolean: true if the BB in any loop 
function BBsInLoop(lives, loop_info)
    in_loop = Dict{Int, Bool}()
    for (j,bb) in lives.basic_blocks
        in_loop[bb.label] = false
    end

    for L in loop_info.loops
        for bbnum in L.members
            in_loop[bbnum] = true
        end
    end
    return in_loop
end

function findMmread(lives, in_loop)
    # Look for such a statement outside any loop: e.g. A = ((top(getfield))(MatrixMarket,:mmread))("b.mtx")
    for (j,bb) in lives.basic_blocks
        if in_loop[bb.label] 
            continue
        end
        for stmt_idx = 1 : length(bb.statements)
            expr = bb.statements[stmt_idx].expr
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
    from_BB       :: Any #LivenessAnalysis.BasicBlock #ISSUE: Not sure why, always complain "LivenessAnalysis" unidentified
    from_stmt_idx :: Int
    to_BB         :: Any #LivenessAnalysis.BasicBlock
    to_stmt_idx   :: Int
end

# One part of a region. Both from/to statements and the other statements between them are included
type RegionInterval
    BB            :: Any #LivenessAnalysis.BasicBlock
    from_stmt_idx :: Int
    to_stmt_idx   :: Int
end

type Region
    mmread_BB       :: Any #LivenessAnalysis.BasicBlock
    mmread_stmt_idx :: Int
    exits           :: Set{ExitEdge}
    intervals       :: Set{RegionInterval}
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
function DFSGrowRegion(mmread_BB, current_BB, start_stmt_idx, lives, in_loop, visited, has_loop, exits, intervals, symbolInfo, loop_info)
    if (DEBUG_LVL >= 2)
        println("\n\nDFSGrowRegion from BB ", current_BB.label, " stmt ", start_stmt_idx, "(", 
            1 <= start_stmt_idx && start_stmt_idx <= length(current_BB.statements) ? 
            current_BB.statements[start_stmt_idx].index : "", ")");
    end
    
    # Disable the requirement of having a loop.
    # TODO: enable it in future
    loop_expected_on_every_path = false
    
    visited[current_BB.label]  = true
    has_loop[current_BB.label] = in_loop[current_BB.label]
    last_stmt_idx = length(current_BB.statements)
    for stmt_idx = start_stmt_idx : last_stmt_idx
        expr = current_BB.statements[stmt_idx].expr
        if expr.head == :return
            if loop_expected_on_every_path && !has_loop[current_BB.label]
                return false
            end
            exit = ExitEdge(current_BB, stmt_idx - 1, current_BB, stmt_idx)
            push!(exits, exit)
            interval = RegionInterval(current_BB, start_stmt_idx, stmt_idx - 1)
            push!(intervals, interval)
            return true
        end
        if expr.head == :throw
            throw("DFSGrowRegion: throw not handled")
            return false
        end
        distributive = checkDistributivity(expr, symbolInfo, true)
        if !distributive
            if loop_expected_on_every_path && !has_loop[current_BB.label]
                return false
            end
            exit = ExitEdge(current_BB, stmt_idx - 1, current_BB, stmt_idx)
            push!(exits, exit)
            interval = RegionInterval(current_BB, start_stmt_idx, stmt_idx - 1)
            push!(intervals, interval)
            return true
        end
    end
    interval = RegionInterval(current_BB, start_stmt_idx, last_stmt_idx) 
    push!(intervals, interval)
    
    # Current BB's statements have been scanned. Now successors
    for succ_BB in current_BB.succs
        if (DEBUG_LVL >= 2)
            println("looking at succ BB ", succ_BB.label, " whose dominators are ", loop_info.dom_dict[succ_BB.label])
        end
        
        if succ_BB.label == -2 # the pseudo exit of the function
            if loop_expected_on_every_path && !has_loop[current_BB.label]
                return false
            end
            exit = ExitEdge(current_BB, last_stmt_idx + 1, succ_BB, 0)
            push!(exits, exit)
            continue
        end            
        
        if !in(mmread_BB.label, loop_info.dom_dict[succ_BB.label])
            # mmread BB does not dominate succ BB
            if loop_expected_on_every_path && !has_loop[current_BB.label]
                return false
            end
            if in_loop[succ_BB.label]
                # we are going to insert reverse reodering on a new block between current and succ BB.
                # If succ BB is in a loop, the new block is also in a loop. We cannot affort that cost 
                return false
            end
            
            exit = ExitEdge(current_BB, last_stmt_idx + 1, succ_BB, 0)
            push!(exits, exit)
            continue
        end
        
        if (visited[succ_BB.label])
            has_loop[current_BB.label] = has_loop[current_BB.label] || has_loop[succ_BB.label]
            continue
        end
        
        success = DFSGrowRegion(mmread_BB, succ_BB, 1, lives, in_loop, visited, has_loop, exits, intervals, symbolInfo, loop_info)
        if !success
            return false
        end
    end
    return true
end

function growRegion(mmread_BB, mmread_stmt_idx, lives, in_loop, symbolInfo, loop_info)
    visited  = Dict{Int, Bool}() # BB has been visited?
    has_loop = Dict{Int, Bool}() # Every path starting from the BB crosses a loop?
    for (j,bb) in lives.basic_blocks 
        visited[bb.label]  = false
        has_loop[bb.label] = false  
    end
    
    # Imagine there are always two invisible statements in any BB, whose statement indices are 0 and last_stmt_idx+1
    # So any real statement always has a statement before and after it. That is why we can do mmread_stmt_idx + 1 here
    # Similary, we can do any stmt_indx -1 as well.
    region = Region(mmread_BB, mmread_stmt_idx, Set{ExitEdge}(), Set{RegionInterval}())
    success = DFSGrowRegion(mmread_BB, mmread_BB, mmread_stmt_idx + 1, lives, in_loop, visited, has_loop, region.exits, region.intervals, symbolInfo, loop_info)

    if (DEBUG_LVL >= 2)
        println("DFSGrowRegion successful?: ", success)
        println("Intervals of region found:")
        for interval in region.intervals
            println("\n\nBB ", interval.BB.label, " stmts ", interval.from_stmt_idx, "(", 
                1 <= interval.from_stmt_idx && interval.from_stmt_idx <= length(interval.BB.statements) ?
                interval.BB.statements[interval.from_stmt_idx].index : "", ") : ", interval.to_stmt_idx, "(",
                1 <= interval.to_stmt_idx && interval.to_stmt_idx <= length(interval.BB.statements) ?
                interval.BB.statements[interval.to_stmt_idx].index : "", ")")
        end
    end
    
    if success
        return region
    else
        return nothing
    end
end

function regionFormationBasedOnMmread(lives, loop_info, symbolInfo)
    in_loop = BBsInLoop(lives, loop_info)
    mmread_BB, mmread_stmt_idx = findMmread(lives, in_loop)
    if mmread_BB == nothing
        return nothing
    end
    return growRegion(mmread_BB, mmread_stmt_idx, lives, in_loop, symbolInfo, loop_info)
end

function findIAs(region::Region, symbolInfo)
    # Build inter-dependent arrays
    IAs = Dict{Any, Set}()
    for interval in region.intervals
        BB = interval.BB
        for stmt_index = interval.from_stmt_idx : interval.to_stmt_idx
            stmt = BB.statements[stmt_index]
            IAs[stmt] = Set{Any}()
            seeds = Set{Any} ()
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

function regionTransformation(funcAST, lives, loop_info, symbolInfo, region, IAs)
    mmread_stmt = region.mmread_BB.statements[ region.mmread_stmt_idx] 
    lhs = mmread_stmt.expr.args[1]

    reordered = Set{Symbol}()
    push!(reordered, lhs)
    
    changed = true
    while changed
        changed = false
        for interval in region.intervals
            BB = interval.BB
            for stmt_idx in interval.from_stmt_idx : interval.to_stmt_idx
                stmt = BB.statements[stmt_idx]
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
    
    # The symbols that should be reordered before the region are the reordered uses live into the region
    reorderedBeforeRegion = intersect(reordered, mmread_stmt.live_out)
    
    dprintln(2, "To be reordered before region: ", reorderedBeforeRegion)

    if isempty(reorderedBeforeRegion)
        return
    end
    
    #TODO: benefit-cost analysis

    # Only those arrays that are updated in the region need to be reverse reordered afterwards
    # Limitation: the region must be linear (except with some loop back edges) in order for the
    # analysis to be correct: if there is any control flow, we might be wrong: some arrays in 
    # the if/else block might never execute, but we assume they do, and propagate their "reordered"
    # info to other arrays. 
    # TODO: minimal instrument to control flow graph to record which arrays are actually updated dynamically
    # if there is control flow. Then insert code to reverse reorder those actually updated.
    updatedInRegion = Set{Symbol}()
    for interval in region.intervals
        BB = interval.BB
        for stmt_idx in interval.from_stmt_idx : interval.to_stmt_idx
            stmt = BB.statements[stmt_idx]
            union!(updatedInRegion, stmt.def)
        end
    end
    dprintln(2, "updatedInRegion: ", updatedInRegion)
      
    # New a vector to hold the new reordering statements R(LiveIn) before the region. 
    new_stmts_before_region = Expr[]

    # Replace the original "lhs=mmread(filename)" with
    #   T   = mmread_reorder(filename) // T is a tuple
    #   lhs = (top(tupleref))(T, 1)
    #   P   = (top(tupleref))(T, 2)
    #   P'  = (top(tupleref))(T, 3)
    T = gensym("T")
    P = gensym("P")
    Pprime = gensym("Pprime")
    
    mmread_stmt.expr.args[1] = T
    mmread_stmt.expr.args[2].args[1].args[3] = QuoteNode(:mmread_reorder) # change mmread to mmread_reorder
    
    stmt = Expr(:(=), lhs, Expr(:call, TopNode(:tupleref), T, 1))
    push!(new_stmts_before_region, stmt)
    
    stmt = Expr(:(=), P, Expr(:call, TopNode(:tupleref), T, 2))
    push!(new_stmts_before_region, stmt)
    
    stmt = Expr(:(=), Pprime, Expr(:call, TopNode(:tupleref), T, 3))
    push!(new_stmts_before_region, stmt)

    # Now reorder other arrays
    for sym in reorderedBeforeRegion
        if sym != lhs
            if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                reorderMatrix(sym, P, Pprime, new_stmts_before_region, false, true, true)
            else
                reorderVector(sym, P, new_stmts_before_region)
            end
        end
    end

    i = 1
    for new_stmt in new_stmts_before_region
        insert!(region.mmread_BB.statements, region.mmread_stmt_idx + i, LivenessAnalysis.TopLevelStatement(0, new_stmt))
        i += 1
    end
                 
    # At each exit edge, we need to reverse reorder those that live out.
    # to recover their original order before they are getting used outside of the region
    reorderedAndUpdated = intersect(reordered, updatedInRegion)
    for exit in region.exits
        if exit.from_BB == exit.to_BB
            BB = exit.from_BB
            # insert statements inside a block. At least one of the
            # from/to stmt_idx is a real statement's index
            last_stmt_idx = length(BB.statements)
            assert((1 <= exit.from_stmt_idx && exit.from_stmt_idx <= last_stmt_idx) ||
                (1 <= exit.to_stmt_idx && exit.to_stmt_idx <= last_stmt_idx))
            if (1 <= exit.from_stmt_idx && exit.from_stmt_idx <= last_stmt_idx)
                stmt = BB.statements[exit.from_stmt_idx]
                reverseReordered = intersect(reorderedAndUpdated, stmt.live_out)
            else 
                stmt = BB.statements[exit.to_stmt_idx]
                reverseReordered = intersect(reorderedAndUpdated, stmt.live_in)
            end
            if isempty(reverseReordered)
                continue
            end
            dprintln(2, "ReverseReorder inside BB ", BB.label, " ", exit.from_stmt_idx, "--> ", exit.to_stmt_idx)
                
            new_stmts = Expr[]
            for sym in reverseReordered
                if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                    reverseReorderMatrix(sym, P, Pprime, new_stmts, false, true, true)
                else
                    reverseReorderVector(sym, P, new_stmts)
                end
            end

            if (1 <= exit.from_stmt_idx && exit.from_stmt_idx <= last_stmt_idx)
                for new_stmt in new_stmts
                    insert!(BB.statements, exit.from_stmt_idx + 1, LivenessAnalysis.TopLevelStatement(0, new_stmt))
                end
            else
                for new_stmt in new_stmts
                    insert!(BB.statements, exit.to_stmt_idx, LivenessAnalysis.TopLevelStatement(0, new_stmt))
                end
            end                
        else
            # insert reverse reordering between two blocks.
            reverseReordered = intersect(reorderedAndUpdated, exit.to_BB.live_in)
            if isempty(reverseReordered)
                continue
            end
            dprintln(2, "ReverseReorder on edge ", exit.from_BB.label, " --> ", exit.to_BB.label)
                
            new_stmts = Expr[]
            for sym in reverseReordered
                if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                    reverseReorderMatrix(sym, P, Pprime, new_stmts, false, true, true)
                else
                    reverseReorderVector(sym, P, new_stmts)
                end
            end
            
            # Now actually change the CFG.
            (new_bb, new_goto_stmt) = LivenessAnalysis.insertBetween(lives, exit.from_BB.label, exit.to_BB.label)
            for new_stmt in new_stmts
                LivenessAnalysis.addStatementToEndOfBlock(lives, new_bb, stmt)
            end
            if new_goto_stmt != nothing
                push!(new_bb.statements, new_goto_stmt)
            end
        end
    end

    if (DEBUG_LVL >= 2)
        println("******** CFG after mmread_reorder: ********")
        show(lives);
    end 
end

function reorderDuringMmread(funcAST, lives, loop_info, symbolInfo)
    if (DEBUG_LVL >= 2)
        println("******** CFG before mmread_reorder: ********")
        show(lives);
    end 
    
    region = regionFormationBasedOnMmread(lives, loop_info, symbolInfo)
    if region == nothing
        return false
    end
    
    IAs = findIAs(region, symbolInfo)
    regionTransformation(funcAST, lives, loop_info, symbolInfo, region, IAs)
    return true
end

function reorder(funcAST, lives, loop_info, symbolInfo)
    assert(funcAST.head == :lambda)
    args = funcAST.args
    assert(length(args) == 3)

    # First, see if there is an mmread() outside all loops.
    # If not, select a sparse matrix from the function AST's arguments. 
    # So far, choose the first sparse matrix in the arguments.
    # TODO: have a heuristic to choose the best candidate? Or explicitly
    # point out the candidate with user annotation?
    # TODO: it really does not matter whether the sparse matrix is from
    # an mmread or from an argument -- from the perspective of region formation,
    # identifying IAs, and region transformation. So should unify the two cases.
    success = reorderDuringMmread(funcAST, lives, loop_info, symbolInfo)
    if (success)
        body_reconstructed = LivenessAnalysis.createFunctionBody(lives)
        funcAST.args[3] = body_reconstructed
        return funcAST
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
        distributive = checkDistributivity(funcAST, symbolInfo, true)
        dprintln(3,"After our type inference, distributive = ", distributive)
          
        if !distributive
            return funcAST
        end
        
        for L in loop_info.loops
            reorderLoop(L, M, lives, symbolInfo)
        end
    end
    
    body_reconstructed   = CompilerTools.LivenessAnalysis.createFunctionBody(lives)
    funcAST.args[3].args = body_reconstructed
    
    funcAST
end

# Try to insert knobs to an expression
function insert_knobs_to_statement(expr:: Expr, lives, loop_info)
    return expr
end

# Analyze the sparse library function calls in the AST of a function, 
# and insert knobs(context-sensitive info)
function insert_knobs(ast, lives, loop_info)	
  if isa(ast, LambdaStaticData)
      ast = uncompressed_ast(ast)
  end
  assert(isa(ast, Expr) && is(ast.head, :lambda))
  body = ast.args[3]
  assert(isa(body, Expr) && is(body.head, :body))
  for i = 1:length(body.args)
	insert_knobs_to_statement(body.args[i], lives, loop_info)    
  end
end

# Try reordering, then insert knobs.
# TODO: combine the two together so that we scan AST only once.
function sparse_analyze(ast, lives, loop_info, symbolInfo)
  dprintln(2, "***************** Sparse analyze *****************")

  ast1 = reorder(ast, lives, loop_info, symbolInfo)
#  ast2 = insert_knobs(ast1, lives, loop_info)
  return ast1
end
