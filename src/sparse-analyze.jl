export CSR_ReorderMatrix, reorderVector, reverseReorderVector

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

# In reordering, we insert some calls to the following 3 functions. So they are executed secretly
function CSR_ReorderMatrix(A::SparseMatrixCSC, newA::SparseMatrixCSC, P::Vector, Pprime::Vector, getPermuation::Bool)
  println("******* Reordering matrix!")
  ccall((:CSR_ReorderMatrix, "../lib/libcsr1.so"), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(newA.colptr), pointer(newA.rowval), pointer(newA.nzval),
               pointer(P), pointer(Pprime), getPermuation)
end

function reorderVector(V::Vector, newV::Vector, P::Vector)
    println("******* Reordering vector!")

   ccall((:reorderVector, "../lib/libcsr1.so"), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         pointer(V), pointer(newV), pointer(P), length(V))
end

function reverseReorderVector(V::Vector, newV::Vector, P::Vector)
   println("******* Reverse reordering vector!")

   ccall((:reorderVectorWithInversePerm, "../lib/libcsr1.so"), Void,
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

# A is the matrix to be reordered. M is the seed matrix we have chosen for doing reorderig 
# analysis. If A and M are the same, we also compute permutation and inverse permutation
# info (P and Pprime). Otherwise, the info has already been computed. That means, this
# function must be called to reorder matrix M first, and then for other matrices
function reorderMatrix(A::Symbol, M::Symbol, P::Symbol, Pprime::Symbol, new_stmts)
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
    getPermuation = (A == M) ? true : false
    stmt = :(CSR_ReorderMatrix($A, $newA, $P, $Pprime, $getPermuation))
    push!(new_stmts, stmt)
    
    # Update the original matrix with the new data. Note: assignment between
    # two arrays actually make them aliased. So there is no copy cost
    stmt = :( $A = $newA )
    push!(new_stmts, stmt)
end

reverseReorderMatrix(sym, M, P, Pprime, landingPad) = 
       reorderMatrix(sym, M, Pprime, P, landingPad)

function reorderVector(V::Symbol, P::Symbol, new_stmts)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the element type. Should make it general in future
    newV = gensym(string(V))
    stmt = :( $newV = Array(Cdouble, size($V, 1))) 
    push!(new_stmts, stmt)
    
    # Do the actual reordering in the C library
    stmt = :(reorderVector($V, $newV, $P))
    push!(new_stmts, stmt)
    
    # Update the original vector with the new data
    stmt = :( $V = $newV )
    push!(new_stmts, stmt)
end

function reverseReorderVector(V::Symbol, P::Symbol, new_stmts)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the element type. Should make it general in future
    newV = gensym(string(V))
    stmt = :( $newV = Array(Cdouble, size($V, 1))) 
    push!(new_stmts, stmt)
    
    # Do the actual reordering in the C library
    stmt = :(reverseReorderVector($V, $newV, $P))
    push!(new_stmts, stmt)
    
    # Update the original vector with the new data
    stmt = :( $V = $newV )
    push!(new_stmts, stmt)
end

function effectedByReordering(S, symbolInfo)
    s1 = Set{Symbol}()
    for sym in S
        typ = symbolInfo[sym]
        if typ <: AbstractArray
            push!(s1, sym)
        end
    end
    s1
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
function reorder(funcAST, L, M, lives, symbolInfo)
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
    # Let "uses/defs(S)" be the set of matrices/vectors used/defined in 
    # statement S in loop L.  
    # Let us say "reorderedUses/Defs" are the matrices/vectors 
    # used/defined in L that got reordered. Then
    #    reorderedUses = union of uses(S) forall statement S 
    #       if uses(S) contains M or any element in reorderedDefs
    #    reorderedDefs = union of defs(S) forall statement S
    #       if uses(S) contains any element in reorderedUses

    bbs = lives.basic_blocks
    
    # In liveness info, we need only matrices and vectors. Build this info
    uses = Dict{Any, Set}()
    defs = Dict{Any, Set}()
    for bbnum in L.members
        for stmt in bbs[bbnum].statements
            uses[stmt] = effectedByReordering(stmt.use, symbolInfo)
            defs[stmt] = effectedByReordering(stmt.def, symbolInfo)
        end
    end
    
    reorderedUses = Set{Symbol}()
    reorderedDefs = Set{Symbol}()

    changed = true
    while changed
        changed = false
        for bbnum in L.members
            for stmt in bbs[bbnum].statements
                if ⊈(uses[stmt], reorderedUses)
                    if in(M, uses[stmt]) || !isempty(intersect(uses[stmt], reorderedDefs))
                        union!(reorderedUses, uses[stmt])
                        changed = true
                    end
                end
                if ⊈(defs[stmt], reorderedDefs)
                    if !isempty(intersect(uses[stmt], reorderedUses))
                        union!(reorderedDefs, defs[stmt])
                        changed = true
                    end
                end
            end
        end
    end
    
    # The symbols that should be reordered before L are the reorderedUses live into L
    headBlock = bbs[L.head]
    reorderedBeforeL = intersect(reorderedUses, headBlock.live_in)
    
    dprintln(2, "To be reordered before L: ", reorderedBeforeL)
    dprintln(2, "To be reordered inside L: ", reorderedDefs)

    if isempty(reorderedBeforeL)
        return
    end
    
    # The symbols that are reordered inside L, due to the symbols reordered before L,
    # are simply the reorderedDefs
    
    # Now we know all the symbols that need/will be reordered
    reordered = union(reorderedBeforeL, reorderedDefs)
           
    #TODO: benefit-cost analysis
      
    # New a vector to hold the new reordering statements R(LiveIn) before the loop. 
    # We do not insert them into the CFG at this moment, since that will change the
    # pred-succ and live info, leading to some subtle errors. We need CFG not changed
    # until all new statements are ready.
    new_stmts_before_L = Expr[]

    # Allocate space to store the permutation and inverse permutation info
    (P, Pprime) = allocateForPermutation(M, new_stmts_before_L)

    # Compute P and Pprime, and reorder M
    reorderMatrix(M, M, P, Pprime, new_stmts_before_L)
    
    # Now reorder other arrays
    for sym in reorderedBeforeL
        if sym != M
            if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                reorderMatrix(sym, M, P, Pprime, new_stmts_before_L)
            else
                reorderVector(sym, P, new_stmts_before_L)
            end
        end
    end
        
    # At the exit of the loop, we need to reverse reorder those that live out of the loop,
    # to recover their original order before they are getting used outside of the loop
    # We remember those in an array. Each element is a tuple 
    # (bb label, succ label, the new statements to be inserted on the edge from bb to succ)
    new_stmts_after_L =  Any[]
    for bbnum in L.members
        bb = bbs[bbnum]
        for succ in bb.succs
            if !in(succ.label, L.members)
                reverseReordered = intersect(reordered, succ.live_in)
                if isempty(reverseReordered)
                    continue
                end
                
                new_stmts = (bbnum, succ.label,  Expr[])
                push!(new_stmts_after_L, new_stmts)
                
                dprintln(2, "ReverseReorder on edge ", bbnum, " --> ", succ.label)
                
                for sym in reverseReordered
                    if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                        reverseReorderMatrix(sym, M, P, Pprime, new_stmts)
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
            reorder(funcAST, L, M, lives, symbolInfo)
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
