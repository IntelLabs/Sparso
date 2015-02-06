# Try to reorder sparse matrices and related (dense vectors) in the given loop
# L is a loop region. M is a sparse matrix live into L. 
# Liveness is the set of matrices or vectors live into/out of L
# Aliases contains pairs of matrices or vectors that alias with each other
# Note: if a matrix/vector is to be reordered, then we must ensure all its 
# definitions reordered in the same way.
# A matrix A is reordered by row and column permuation, i.e. P'AP, where P is
# a permuation matrix, and P' is its conjugate transpose.
# A (column) vector v is reordered by row permutation, i.e. P'v
function reorder(funcAST, L, M::Symbol, liveness, uniqSet)
    # If we reorder M inside L, consequently, some other matrices or vectors
    # need to be reordered. For example, for the following code
    #       N = M * v 
    #       K = N + speye(size(N, 1))
    # v has to be reordered too, to make M * v meaningful. As the result of 
    # M * v, N got reordered as well. Consequently, the matrix due to speye()
    # has to be reordered, and K as well. This is a propagation.
    # Assumption: we know the semantics of all the operators/functions.
    # Let "uses/defs(S)" be the set of matrices/vectors used/defined in the 
    # RHS/LHS of a statement S in loop L.  
    # Let us say "reorderedUses/Defs" are the 
    # matrices/vectors used/defined in L that got reordered. Then
    #    reordedUses = union of uses(S) forall statement S 
    #       if uses(S) contains M or any element in reorderedDefs
    #    reordedDefs = union of defs(S) forall statement S
    #       if uses(S) contains any element in reorderedUses
  
    reordedUses = Set{Symbol}()
    reordedDefs = Set{Symbol}()

    bbs = lives.basic_blocks
    changed = true
    while changed
        changed = false
        for bbnum in L
            for stmt in bbs[bbnum].statements
                if ⊈(stmt.use, reordedUses)
                    if in(M, stmt.use) or !isempty(intersect(stmt.use, reordedDefs))
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
    
    # Permuted_Before_L is the reorderedUses live into L
    # Permuted_Inside_L is just reorderedDefs
Permuted_Before_L =  { matrix or vector live into L and in an expression with M transitively}
Permuted_Inside_L =  { matrix or vector defined in L depending on Permuted_Before_L  transitively}
Permuted = Permuted_Before_L  |  Permuted_Inside_L
Liveout = {matrix or vector live out of L} //Todd
ReversePermuted = Permuted & Liveout
E     = {expression in L containing a matrix or vector in Permuted }
Insert:
             Benefit = benefit of computing E with Permuted 
             Cost = cost of permuting Permuted_Before_L and reverse permuting with ReversePermuted
             If (benefit > cost) {
                                                      Before loop: insert “permute Permuted_Before_L" "ccall(reordering, librarypath, signature)
                                                      after loop: insert “reverse permute ReversePermuted”
             }

Expressions are limited to the following form and their combinations: 
             A*v, A+-A, A*A, v*v, v+-v, n*v, n*n, n+-n
Where A: matrix, v: vector, n: number 

end

# Try to insert knobs to an expression
insert_knobs_to_expr(expr:: Expr, invariants) =
	expr.head == :(call) ? insert_knobs_to_call(expr, invariants) : 
	is_assignment(expr.head) ? insert_knobs_to_assignment(expr, invariants) :
	expr.head == :(comparison) ? insert_knobs_to_comparison(expr, invariants) :
	expr.head == :(block) ? insert_knobs_to_block(expr, invariants) 
ex.head == :(.) ? tqvar(ex) :
		ex.head == :(call) ? tcall(ex) :
		ex.head == :(comparison) ? tcomparison(ex) :
		ex.head == :(ref) ? tref(ex) :
		ex.head == :(=) ? tassign(ex) :
		ex.head == :(block) ? tblock(ex) :
		is_empty_tuple(ex) ? TEmpty() :
		is_opassign(ex.head) ? topassign(ex) :

	

# Analyze the sparse library function calls in the AST of a function, 
# and insert knobs(context-sensitive info)
function insert_knobs(ast, loop, invariants)	
  if isa(ast, LambdaStaticData)
      ast = uncompressed_ast(ast)
  end
  assert(isa(ast, Expr) && is(ast.head, :lambda))
  body = ast.args[3]
  assert(isa(body, Expr) && is(body.head, :body))
  for i = 1:length(body.args)
	insert_knobs_to_statement(body.args[i], invariants)    
  end
end

function sparse_analyze(ast, lives, loop_info)
  dprintln(2, "sparse_analyze: ", ast, " ", lives, " ", loop_info)
  ast1 = reorder(ast, lives, loop_info, invariants)
  ast2 = insert_knobs(ast1, lives, loop_info, invariants)
  return new_ast2
end
