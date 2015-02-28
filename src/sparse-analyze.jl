export CSR_ReorderMatrix, reorderVector, reverseReorderVector

# In reordering, we insert some calls to the following 3 functions. So they are executed secretly
function CSR_ReorderMatrix(A::SparseMatrixCSC, newA::SparseMatrixCSC, P::Vector, Pprime::Vector, getPermuation::Bool)
  ccall((:CSR_ReorderMatrix, "../lib/libcsr.so"), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Cbool),
               A.m, A.n, A.colptr, A.rowval, A.nzval, 
               newA.colptr, newA.rowval, newA.nzval,
               P, Pprime, getPermuation)
end

function reorderVector(V::Vector, newV::Vector, P::Vector)
   ccall((:reorderVector, "../lib/libcsr.so"), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         V, newV, P, length(V))
end

function reverseReorderVector(V::Vector, newV::Vector, Pprime::Vector)
   ccall((:reorderVectorWithInversePerm, "../lib/libcsr.so"), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
         V, newV, Pprime, length(V))
end

function allocateForPermutation(M::Symbol, lives, block)
    P = gensym("P")
    stmt = :($P = Array(Cint, size($M, 2)))
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    Pprime = gensym("Pprime")
    stmt = :($Pprime = Array(Cint, size($M, 2)))
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    (P, Pprime)
end

# A is the matrix to be reordered. M is the seed matrix we have chosen for doing reorderig 
# analysis. If A and M are the same, we also compute permutation and inverse permutation
# info (P and Pprime). Otherwise, the info has already been computed. That means, this
# function must be called to reorder matrix M first, and then for other matrices
function reorderMatrix(A::Symbol, M::Symbol, P::Symbol, Pprime::Symbol, lives, block)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the sparse matrix format and element type. Should make it
    # general in future
    newA = gensym(string(A))
    stmt = :(  
        $newA = SparseMatrixCSC($A.m, $A.n, 
                  Array(Cint, size($A.colptr, 1)), 
                  Array(Cint, size($A.rowval, 1)), 
                  Array(Cdouble, size($A.nzval, 1))) 
    )
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    # Do the actual reordering in the C library
    getPermuation = (A == M) ? true : false
    stmt = :(CSR_ReorderMatrix($A, $newA, $P, $Pprime, $getPermuation))
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    # Update the original matrix with the new data. Note: assignment between
    # two arrays actually make them aliased. So there is no copy cost
    stmt = :( $A = $newA )
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
end

                        
reverseReorderMatrix(sym, M, P, Pprime, lives, landingPad) = 
       reorderMatrix(sym, M, Pprime, P, lives, landingPad)

function reorderVector(V::Symbol, P::Symbol, lives, block)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the element type. Should make it general in future
    newV = gensym(string(V))
    stmt = :( $newV = Array(Cdouble, size($V, 1))) 
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    # Do the actual reordering in the C library
    stmt = :(reorderVector($V, $newV, $P))
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    # Update the original vector with the new data
    stmt = :( $V = $newV )
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
end

function reverseReorderVector(V::Symbol, Pprime::Symbol, lives, block)
    # Allocate space that stores the reordering result in Julia
    # TODO: Here we hard code the element type. Should make it general in future
    newV = gensym(string(V))
    stmt = :( $newV = Array(Cdouble, size($V, 1))) 
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    # Do the actual reordering in the C library
    stmt = :(reorderVectorWithInversePerm($V, $newV, $Pprime))
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
    
    # Update the original vector with the new data
    stmt = :( $V = $newV )
    LivenessAnalysis.addStatementToEndOfBlock(lives, block, stmt)
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
        
    # Insert R(LiveIn) before the head. If the head has more than,
    # one predecessor, insert an empty block before the head.
    landingPad = nothing
    if length(headBlock.preds) == 1
        landingPad = first(headBlock.preds)
    else
        landingPad = LivenessAnalysis.insertBefore(lives, L.head)
    end

    # Allocate space to store the permutation and inverse permutation info
    (P, Pprime) = allocateForPermutation(M, lives, landingPad)

    # Compute P and Pprime, and reorder M
    reorderMatrix(M, M, P, Pprime, lives, landingPad)
    
    # Now reorder other arrays
    for sym in reorderedBeforeL
        if sym != M
            if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                reorderMatrix(sym, M, P, Pprime, lives, landingPad)
            else
                reorderVector(sym, P, lives, landingPad)
            end
        end
    end
    
    if(DEBUG_LVL >= 2)
        println("******** Landing pad before L ********")
        show(landingPad);
    end 
    
    # At the exit of the loop, we need to reverse reorder those that live out of the loop,
    # to recover their original order before they are getting used outside of the loop
    for bbnum in L.members
        bb = bbs[bbnum]
        for succ in bb.succs
            if !in(succ.label, L.members)
                reverseReordered = intersect(reordered, succ.live_in)
                if isempty(reverseReordered)
                    continue
                end
                landingPad = nothing
                if length(succ.preds) == 1
                    landingPad = succ
                else
                    landingPad = LivenessAnalysis.insertBetween(lives, bbnum, succ.label)
                end
                
                for sym in reverseReordered
                    if typeOfNode(sym, symbolInfo) <: AbstractMatrix
                        reverseReorderMatrix(sym, M, P, Pprime, lives, landingPad)
                    else
                        reverseReorderVector(sym, Pprime, lives, landingPad)
                    end
                end
                
                if(DEBUG_LVL >= 2)
                    println("******** Landing pad after L for bb ", bbnum, " ********")
                    show(landingPad);
                end 
            end
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
