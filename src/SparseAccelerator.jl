module SparseAccelerator

export @acc

# Sparse ACcelerator entry. 
macro acc(ast_node)
  foreach callsite cs in ast_node
    cs_ast = code_typed(cs, sig)
	loops = existing_liveness_code(cs_ast)
	invariants = existing_alias_analysis(cs_ast, loops)
	#TODO: this will be changed to a Julia call, which seems more convenient than C in manipulating AST
	new_cs_name = ccall((:analyze, libpsl_analyzer), Ptr{Void}, (Ptr{Void}, Ptr{Void}, Ptr{Void}), ast_node, loops, invariants) 
	new_expr = new_expr + new_cs_name
  end
  new_expr
end

