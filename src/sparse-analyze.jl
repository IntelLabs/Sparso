export CSR_ReorderMatrix, reorderVector, reverseReorderVector

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

# In reordering, we insert some calls to the following 3 functions. So they are executed secretly
# Reorder sparse matrix A and store the result in newA. A itself is not changed.
function CSR_ReorderMatrix(A::SparseMatrixCSC, newA::SparseMatrixCSC, P::Vector, Pprime::Vector, getPermuation::Bool, oneBasedInput::Bool, oneBasedOutput::Bool)
  ccall((:CSR_ReorderMatrix, "../lib/libcsr.so"), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool, Bool, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(newA.colptr), pointer(newA.rowval), pointer(newA.nzval),
               pointer(P), pointer(Pprime), getPermuation, oneBasedInput, oneBasedOutput)
end

function CSR_Bandwidth(A::SparseMatrixCSC)
   A2 = ccall((:CSR_Create, "../lib/libcsr.so"), Ptr{Void},
         (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint),
         A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval), 1)
   bw = ccall((:CSR_GetBandwidth, "../lib/libcsr.so"), Cint,
         (Ptr{Void},),
         A2)
   ccall((:CSR_Destroy, "../lib/libcsr.so"), Void,
         (Ptr{Void},),
         A2)
   bw
end

function reorderVector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVector, "../lib/libcsr.so"), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(newV), pointer(P), length(V))
end

function reverseReorderVector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVectorWithInversePerm, "../lib/libcsr.so"), Void,
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
          ast.args[1] = LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!))
          return ast
        end
      elseif length(ast.args) == 6
        alpha = ast.args[2]
        A = ast.args[3]
        x = ast.args[4]
        beta = ast.args[5]
        y = ast.args[6]
        if A.typ <: SparseMatrixCSC && x.typ <: Vector
          ast.args[1] = LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!))
          return ast
        end
      end
    elseif ast.args[1] == :(*)
      dprintln(3,"optimize_calls found *")
      A = ast.args[2]
      x = ast.args[3]
      if A.typ <: SparseMatrixCSC && x.typ <: Vector
        dprintln(3,"optimize_calls converting to SpMV")
        ast.args[1] = LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV))
        return ast
      end
    end
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
            # Replace calls to A_mul_B! with calls to SpMV!
            bbs[bbnum].statements[stmt_index].expr = AstWalker.get_one(AstWalker.AstWalk(bbs[bbnum].statements[stmt_index].expr, optimize_calls, nothing))
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
    (new_bb, new_goto_stmt) = LivenessAnalysis.insertBefore(lives, L.head, true, L.back_edge)
    for stmt in new_stmts_before_L
        LivenessAnalysis.addStatementToEndOfBlock(lives, new_bb, stmt)
    end
    if new_goto_stmt != nothing
      push!(new_bb.statements, new_goto_stmt)
    end
    
    for (pred, succ, new_stmts) in new_stmts_after_L
        (new_bb, new_goto_stmt) = LivenessAnalysis.insertBetween(lives, pred, succ)
        for stmt in new_stmts
            LivenessAnalysis.addStatementToEndOfBlock(lives, new_bb, stmt)
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

function reorder(funcAST, lives, loop_info, symbolInfo)
    assert(funcAST.head == :lambda)
    args = funcAST.args
    assert(length(args) == 3)
    local param = args[1]

    # Select a sparse matrix from the function AST's arguments. 
    # So far, choose the first sparse matrix in the arguments.
    # TODO: have a heuristic to choose the best candidate?
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
        for L in loop_info.loops
            reorderLoop(L, M, lives, symbolInfo)
        end
    end
    
    body_reconstructed = LivenessAnalysis.createFunctionBody(lives)
    funcAST.args[3] = body_reconstructed
    
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
