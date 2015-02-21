function RCM{Tv,Ti<:Integer}(M::SparseMatrixCSC{Tv, Ti})
    #call to Jongsoo's library. Return a permutation matrix P
    #P = ccall((:RCM, pcl libary path), (Int, Int, Vector{Ti}, Vector{Ti}, Vector{Tv}), M.m, M.n, M.colptr, M.rowval, m.nzval)
    P = eye(size(M, 1))
    return P
end


# A matrix A is reordered by row and column permuation, i.e. P'AP, where P is
# a permuation matrix, and P' is its conjugate transpose.
# A (column) vector v is reordered by row permutation, i.e. P'v
# TODO: insert a new statement into the right place of AST and a basic block

function reorder(S::Symbol, P::AbstractMatrix)
    typ = get_sym_type(S)
    if typ :< AbstractSparseMatrix
        return :(M = P' * M * P)
    elseif typ :< AbstractVector
        return :(V = P' * V)
    else
        assert(false)
    end 
end

function reverseReorder(S::Symbol, P::AbstractMatrix)
    typ = get_sym_type(S)
    if typ :< AbstractSparseMatrix
        return :(M = P * M * P')
    elseif typ :< AbstractVector
        return :(V = P * V)
    else
        assert(false)
    end 
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
    
function reorder(funcAST, L, M::AbstractSparseMatrix, lives)
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
    #    reordedUses = union of uses(S) forall statement S 
    #       if uses(S) contains M or any element in reorderedDefs
    #    reordedDefs = union of defs(S) forall statement S
    #       if uses(S) contains any element in reorderedUses
        
    reordedUses = Set{Symbol}()
    reordedDefs = Set{Symbol}()

    changed = true
    while changed
        changed = false
        for bbnum in L.members
            for stmt in bbs[bbnum].statements
                if ⊈(stmt.use, reordedUses)
                    if in(M.name, stmt.use) || !isempty(intersect(stmt.use, reordedDefs))
                        union!(reordedUses, stmt.use)
                        changed = true
                    end
                end
                if ⊈(stmt.def, reorderedDefs)
                    if !isempty(intersect(stmt.use, reordedUses))
                        union!(reordedDefs, stmt.def)
                        changed = true
                    end
                end
            end
        end
    end
    
    # The symbols that should be reordered before L are the reorderedUses live into L
    reorderedBeforeL = intersect(reorderedUses, L.head.live_in)
    
    # The symbols that are reordered inside L, due to the symbols reordered before L,
    # are simply the reorderedDefs
    
    # Now we know all the symbols that need/will be reordered
    reordered = union(reorderedBeforeL, reorderedDefs)
           
    #TODO: benefit-cost analysis
    
    # Get the permutation matrix
    P = RCM(M)
    
    # Insert R(LiveIn) before the head
    for sym in reorderedBeforeL
        reorder(sym, P)
    end
    
    # At the exit of the loop, we need to reverse reorder those that live out of the loop,
    # to recover their original order before they are getting used outside of the loop
    # TODO: change this after Todd add empty basic blocks
    for bbnum in L.exits
        bb = bbs[bbnum]
        reverseReordered = intersect(reordered, bb.live_out)
        
        #insert transformation here
        for sym in reverseReordered
            reverseReorder(sym, P)
        end
    end
end

function reorder(funcAST, lives, loop_info)
    assert(funcAST.head == :lambda)
    args = funcAST.args
    assert(length(args) == 3)
    local param = args[1]

    # Select a sparse matrix from the function AST's arguments. 
    # So far, choose the first sparse matrix in the arguments.
    # TODO: have a heuristic to choose the best candidate?
    found = false
    for i = 1:length(param)
        if typeof(param[i]) <: AbstractSparseMatrix
            found = true
            M = param[i]
            break
        end
    end    
    if found 
        for L in loop_info
            reorder(funcAST, L, M, lives)
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
function sparse_analyze(ast, lives, loop_info)
  dprintln(2, "sparse_analyze: ", ast, " ", lives, " ", loop_info)
  ast1 = reorder(ast, lives, loop_info)
#  ast2 = insert_knobs(ast1, lives, loop_info)
  return ast1
end
