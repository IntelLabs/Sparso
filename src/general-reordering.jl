#=
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

# Reordering analysis: Feasibility, and where and what to reorder.

# ISSUE: We assume that there are NO aliases between any two arrays. We should
# use both static alias analysis and dynamic alias check to get rid of the 
# assumption.

@doc """
A map from a permutation vector to its inverse.
"""
const inverse_color_map = Dict{Int, Int}(
    ROW_PERM     => ROW_INV_PERM,
    ROW_INV_PERM => ROW_PERM,
    COL_PERM     => COL_INV_PERM,
    COL_INV_PERM => COL_PERM
)

@doc """
A vertex in the inter-dependence graph. Each vertex represents an array's
row or column permutation vector. Color indicates the vector.
"""
type InterDependenceGraphVertex
    array      :: Symexpr
    row_perm   :: Bool
    color      :: Int
    neighbours :: Set{Tuple{InterDependenceGraphVertex, Bool}}
    
    InterDependenceGraphVertex(_array, _row_perm) = new(_array, _row_perm, NOT_PERM_YET,
                               Set{Tuple{InterDependenceGraphVertex, Bool}}())
end

@doc """
Inter-dependence graph. It separately represents rows and columns of arrays. It
will be colored starting from the seed.
"""
type InterDependenceGraph
    rows             :: Dict{Symexpr, InterDependenceGraphVertex}
    columns          :: Dict{Symexpr, InterDependenceGraphVertex}
    seed             :: Sym
    perm_restriction :: Int
    
    InterDependenceGraph(_seed) = new(
        Dict{Symexpr, InterDependenceGraphVertex}(),
        Dict{Symexpr, InterDependenceGraphVertex}(),
        _seed, PERM_VECTORS_ARE_INDEPENDENT)
end
 
@doc """
The CallSites' extra field for reordering.
"""
type ReorderingExtra
    seed                     :: Sym
    decider_ast              :: Expr
    matrix2mknob             :: Dict{Symexpr, Symbol}
    decider_bb               :: Any # BasicBlock
    decider_stmt_idx         :: Int
    bb                       :: Any # BasicBlock
    stmt_idx                 :: StatementIndex
    prev_stmt_idx            :: StatementIndex
    inter_dependence_graph   :: InterDependenceGraph

    # Scratch variables.
    live_in_before_prev_expr :: Set{Sym}
    live_in_before_expr      :: Set{Sym}

    ReorderingExtra(_seed, _decider_ast, _matrix2mknob) = new(_seed, _decider_ast,
             _matrix2mknob, nothing, 0, nothing, 0, 0, InterDependenceGraph(_seed),
             Set{Sym}(), Set{Sym}())
end

@doc """
Argument types of setindex!.

A matrix has the following access patterns: [:, j], [i, :], [:, :], and [i, j] 
A vector has the following access patterns: [i], and [:]. 
"""
const setindex!_arg_types1 = (AbstractMatrix, Any, GlobalRef(Main, :(:)), Number)
const setindex!_arg_types2 = (AbstractMatrix, Any, Number, GlobalRef(Main, :(:)))
const setindex!_arg_types3 = (AbstractMatrix, Any, GlobalRef(Main, :(:)), GlobalRef(Main, :(:)))
const setindex!_arg_types4 = (AbstractMatrix, Any, Number, Number)
const setindex!_arg_types5 = (AbstractVector, Any, Number)
const setindex!_arg_types6 = (AbstractVector, Any, GlobalRef(Main, :(:)))

@doc """
For setindex!(A, x, I), map from the setindex! skeleton to the relation
between A and  x. If a dimension of A should not be permuted, then map to
NEVER_PERM.
"""
const setindex!2relations = Dict(
    setindex!_arg_types1 => (ROW_ROW,    NEVER_PERM   ),
    setindex!_arg_types2 => (NEVER_PERM, COLUMN_ROW   ),
    setindex!_arg_types3 => (ROW_ROW,    COLUMN_COLUMN),
    setindex!_arg_types4 => (NEVER_PERM, NEVER_PERM   ),
    setindex!_arg_types5 => (NEVER_PERM,              ),
    setindex!_arg_types6 => (ROW_ROW,                 ),
)

@doc """
Argument types of getindex.

A matrix has the following access patterns: [:, j], [i, :], [:, :], and [i, j] 
A vector has the following access patterns: [i], and [:]. 
"""
const getindex_arg_types1  = (AbstractMatrix, GlobalRef(Main, :(:)), Number)
const getindex_arg_types2  = (AbstractMatrix, Number, GlobalRef(Main, :(:)))
const getindex_arg_types3  = (AbstractMatrix, GlobalRef(Main, :(:)), GlobalRef(Main, :(:)))
const getindex_arg_types4  = (AbstractMatrix, Number, Number)
const getindex_arg_types5  = (AbstractVector, Number)
const getindex_arg_types6  = (AbstractVector, GlobalRef(Main, :(:)))

@doc """
For getindex(A, I), map from the getindex skeleton to the relation
between A and getindex(A, I). If a dimension of A should not be permuted, 
then map to NEVER_PERM.
"""
const getindex2relations = Dict(
    getindex_arg_types1  => (ROW_ROW,    NEVER_PERM   ),
    getindex_arg_types2  => (NEVER_PERM, COLUMN_ROW   ),
    getindex_arg_types3  => (ROW_ROW,    COLUMN_COLUMN),
    getindex_arg_types4  => (NEVER_PERM, NEVER_PERM   ),
    getindex_arg_types5  => (NEVER_PERM,              ),
    getindex_arg_types6  => (ROW_ROW,                 )
)

@doc """"
Build an edge between the two vertices for the inter-dependence graph.
Inverse indicates if their permutation vectors are inverse to each other.
"""
function build_edge(
    vertex1 :: InterDependenceGraphVertex,
    vertex2 :: InterDependenceGraphVertex,
    inverse :: Bool
)
    push!(vertex1.neighbours, (vertex2, inverse))
    push!(vertex2.neighbours, (vertex1, inverse))
end

@doc """
Build a vertex for the array's row or column permutation vector in the
inter-dependence graph.
"""
function build_vertex(
    array    :: Symexpr,
    row_perm :: Bool,
    graph    :: InterDependenceGraph
)
    if row_perm
        if !haskey(graph.rows, array)
            vertex = InterDependenceGraphVertex(array, row_perm)
            graph.rows[array] = vertex
            return vertex
        end
        return graph.rows[array]
    else
        if !haskey(graph.columns, array)
            vertex = InterDependenceGraphVertex(array, row_perm)
            graph.columns[array] = vertex
            return vertex
        end
        return graph.columns[array]
    end
end

@doc """
Build two vertices and and an edge between them in the inter-dependence graph.
"""
function build_vertices_and_edge(
    array1   :: Any,
    array2   :: Any,
    relation :: Int,
    graph    :: InterDependenceGraph
)
    if relation == ROW_ROW
        vertex1 = build_vertex(array1, true, graph)
        vertex2 = build_vertex(array2, true, graph)
        build_edge(vertex1, vertex2, false)
    elseif relation == COLUMN_COLUMN
        vertex1 = build_vertex(array1, false, graph)
        vertex2 = build_vertex(array2, false, graph)
        build_edge(vertex1, vertex2, false)
    elseif relation == COLUMN_ROW_INVERSE
        vertex1 = build_vertex(array1, false, graph)
        vertex2 = build_vertex(array2, true, graph)
        build_edge(vertex1, vertex2, true)
    else
        assert(relation == COLUMN_ROW)
        vertex1 = build_vertex(array1, false, graph)
        vertex2 = build_vertex(array2, true, graph)
        build_edge(vertex1, vertex2, false)
    end
end

@doc """
Build inter-dependence graph with the inter-dependence information drawn from
the AST of a function call to setindex!.
"""
function build_inter_dependence_graph_for_setindex!(
    ast       :: Any,
    arg_types :: Tuple,
    args      :: Vector,
    relations :: Tuple,
    graph     :: InterDependenceGraph,
)
    assert(arg_types[1] <: AbstractMatrix || arg_types[1] <: AbstractVector)
      
    A = get_symexpr(args[1])

    # Build dependence between A and the AST, which represents the result of
    # setindex!(). Since setindex!(A, x, I) returns the updated A, the relations
    # between them are simply ROW_ROW and COLUMN_COLUMN.
    build_vertices_and_edge(A, ast, ROW_ROW, graph)
    if arg_types[1] <: AbstractMatrix
        build_vertices_and_edge(A, ast, COLUMN_COLUMN, graph)
    end

    # Now build dependence between A and x.
    x = get_symexpr(args[2])
    for i in 1 : length(relations)
        relation = relations[i]
        if relation == NEVER_PERM
            # This is not actually a relation between A and x. Instead, it says
            # the first/second dimension of A (A's rows/columns) cannot be permuted.
            if i == 1
                graph.rows[A].color    = NEVER_PERM
            else 
                graph.columns[A].color = NEVER_PERM
            end
        else
            assert(relation == ROW_ROW || relation == COLUMN_ROW || relation == COLUMN_COLUMN)
            build_vertices_and_edge(A, x, relation, graph)
        end
    end
end

@doc """
Build inter-dependence graph with the inter-dependence information drawn from
the AST of a function call to getindex.
"""
function build_inter_dependence_graph_for_getindex(
    ast       :: Any,
    arg_types :: Tuple,
    args      :: Vector,
    relations :: Tuple,
    graph     :: InterDependenceGraph,
)
    assert(arg_types[1] <: AbstractMatrix || arg_types[1] <: AbstractVector)
      
    A = get_symexpr(args[1])

    # Build dependence between A and the AST, which represents the result of
    # getindex().
    for i in 1 : length(relations)
        relation = relations[i]
        if relation == NEVER_PERM
            # This is not actually a relation between A and x. Instead, it says
            # the first/second dimension of A (A's rows/columns) cannot be permuted.
            if i == 1
                graph.rows[A].color    = NEVER_PERM
            else 
                graph.columns[A].color = NEVER_PERM
            end
        else
            assert(relation == ROW_ROW || relation == COLUMN_ROW || relation == COLUMN_COLUMN)
            build_vertices_and_edge(A, ast, relation, graph)
        end
    end
end

@doc """
Build inter-dependence graph with the inter-dependence information drawn from
the AST of a function call, if the call is setindex! or getindex.
"""
function build_inter_dependence_graph_for_set_get_index(
    ast             :: Any,
    module_name     :: AbstractString,
    function_name   :: AbstractString,
    arg_types       :: Tuple,
    args            :: Vector,
    graph           :: InterDependenceGraph,
    handle_getindex :: Bool
)
    if module_name != "Main" || (function_name != "setindex!" && function_name != "getindex") 
        return false
    end

    index2relations = handle_getindex ? getindex2relations : setindex!2relations 
    for (skeleton, relations) in index2relations
        # We cannot use match_skeletons(arg_types, skeleton), since arg_types
        # contains a GlobalRef type, not a GlobalRef(Main, :(:)) (which is not a
        # type, but an instantiation of the type).
        matched = true
        if length(args) == length(skeleton)
            for i in 1 : length(skeleton)
                if skeleton[i] == GlobalRef(Main, :(:))
                    if args[i] != GlobalRef(Main, :(:))
                        matched = false
                        break
                    end
                else
                    if !(arg_types[i] <: skeleton[i])
                        matched = false
                        break
                    end
                end
            end
        else
            matched = false        
        end

        if matched
            if handle_getindex 
                build_inter_dependence_graph_for_getindex(ast, arg_types, args, relations, graph)
            else
                build_inter_dependence_graph_for_setindex!(ast, arg_types, args, relations, graph)
            end
            return true
        end
    end
    return false
end

@doc """
Build inter-dependence graph with the inter-dependence information drawn from
the AST of a function call, if the call is setindex! or getindex.
"""
function build_inter_dependence_graph_for_set_get_index(
    ast             :: Any,
    module_name     :: AbstractString,
    function_name   :: AbstractString,
    arg_types       :: Tuple,
    args            :: Vector,
    graph           :: InterDependenceGraph,
    handle_getindex :: Bool
)
    if module_name != "Main" || (function_name != "setindex!" && function_name != "getindex") 
        return false
    end

    index2relations = handle_getindex ? getindex2relations : setindex!2relations 
    for (skeleton, relations) in index2relations
        # We cannot use match_skeletons(arg_types, skeleton), since arg_types
        # contains a GlobalRef type, not a GlobalRef(Main, :(:)) (which is not a
        # type, but an instantiation of the type).
        matched = true
        if length(args) == length(skeleton)
            for i in 1 : length(skeleton)
                if skeleton[i] == GlobalRef(Main, :(:))
                    if args[i] != GlobalRef(Main, :(:))
                        matched = false
                        break
                    end
                else
                    if !(arg_types[i] <: skeleton[i])
                        matched = false
                        break
                    end
                end
            end
        else
            matched = false        
        end

        if matched
            if handle_getindex 
                build_inter_dependence_graph_for_getindex(ast, arg_types, args, relations, graph)
            else
                build_inter_dependence_graph_for_setindex!(ast, arg_types, args, relations, graph)
            end
            return true
        end
    end
    return false
end

@doc """
Get the symbol or expr representing the array. If array_index is something like 
"0.1", it means the return result of the AST is a tuple, and element 1 of the tuple
is needed; in this case, we return a manually made expression 
Expr(:TupleElement, ast, 1).
"""
function get_array(
    array_index :: Union{Int, AbstractString},
    ast         :: Expr,
    args        :: Vector,
)
    if typeof(array_index) <: Int
        a = (array_index == 0) ? ast :
             (typeof(args[array_index]) == SymbolNode ?
              args[array_index].name : args[array_index])
        return a
    else   
        assert(array_index[1] == '0' && array_index[2] == '.')
        index = parse(Int, array_index[3:end])
        return Expr(:TupleElement, ast, index)
    end
end

@doc """
Build inter-dependence graph with the inter-dependence information drawn from
the AST of a function call.
"""
function build_inter_dependence_graph_for_call(
    ast           :: Any,
    module_name   :: AbstractString,
    function_name :: AbstractString,
    arg_types     :: Tuple,
    args          :: Vector,
    call_sites    :: CallSites
)
    # Recursively handle each argument 
    for arg in args
        build_inter_dependence_graph(arg, call_sites)
    end

    if function_name == "tuple" || function_name == "apply_type"
        # TODO: for each element of a tuple, build an interdependence between
        # it and an Expr(TupleElement, result ast, index of the element)
        return nothing
    end

    # Now handle the call itself
    all_numbers, some_arrays = numbers_or_arrays(ast.typ, arg_types)
    if all_numbers || !some_arrays
        # The function call's result and arguments are all numbers, or 
        # some are non-numbers (like Range{UInt64}) but not regular arrays. 
        # No arrays to care about.
        return nothing
    end

    if function_name == ""
        throw(UnresolvedFunction(ast.head, args[1]))
    end

    graph = call_sites.extra.inter_dependence_graph
    
    # Special handling to setindex! and getindex: so far, we do not generate code
    # for individual array element access. Thus we need to indicate that the
    # index of the array element is not subject to reordering (Set the rows 
    # or columns as NEVER_PERMUTE).
    # Try setindex!
    if build_inter_dependence_graph_for_set_get_index(ast, module_name, function_name, arg_types, args, graph, false)
        return nothing
    end

    # Try getindex
    if build_inter_dependence_graph_for_set_get_index(ast, module_name, function_name, arg_types, args, graph, true)
        return nothing
    end

    fd = look_for_function_description(module_name, function_name, arg_types)
    if fd == nothing
        throw(UndescribedFunction(module_name, function_name, arg_types, ast))
    end
    if !fd.distributive
        throw(NonDistributiveFunction(module_name, function_name, arg_types))
    end

    for (array_index1, array_index2, relation) in fd.IA
        array1 = get_array(array_index1, ast, args)
        array2 = get_array(array_index2, ast, args)
        build_vertices_and_edge(array1, array2, relation, graph)
    end

    return nothing
end
    
@doc """
Build inter-dependence graph with the inter-dependence information drawn from
the AST (of a statement or part of a statement).
"""
function build_inter_dependence_graph(
    ast        :: Any,
    call_sites :: CallSites
)
    if ast == call_sites.extra.decider_ast
        assert(call_sites.extra.decider_bb == nothing)
        assert(call_sites.extra.decider_stmt_idx == 0)
        call_sites.extra.decider_bb       = call_sites.extra.bb
        call_sites.extra.decider_stmt_idx = call_sites.extra.stmt_idx
    end


    local asttyp = typeof(ast)
    if asttyp <: Tuple || asttyp <: Array
        for i = 1:length(ast)
            build_inter_dependence_graph(ast[i], call_sites)
        end
        return nothing
    elseif asttyp == Expr
        symbol_info = call_sites.symbol_info
        head        = ast.head
        args        = ast.args
        if head == :return
            return build_inter_dependence_graph(args, call_sites)
        elseif head == :gotoifnot
            if_clause = args[1]
            return build_inter_dependence_graph(if_clause, call_sites)
        elseif head == :line
            return nothing
        elseif head == :call || head == :call1
            if type_of_ast_node(args[end], symbol_info) == FunctionKnob
                # Ignore the fknob added to the call.
                num_args = length(args) - 2
            else
                num_args = length(args) - 1
            end
            module_name, function_name = resolve_call_names(args)
            arg_types                  = ntuple(i-> type_of_ast_node(args[i+1],
                                                symbol_info), num_args)
            arguments                  = args[2 : num_args + 1]
            return build_inter_dependence_graph_for_call(ast, module_name,
                                function_name, arg_types, arguments, call_sites)
        elseif head == :(=)
            prev_stmt_idx = call_sites.extra.prev_stmt_idx
            if prev_stmt_idx > 0
                # Test if the previous and the current expr compose of a pattern
                # like this:
                #   GenSym(1) = (top(indexed_next))(GenSym(0),x,#s73)    
                #   L = (top(getfield))(GenSym(1),1)
                # where GenSym(0) is a tuple of SparseMatrixCSC.
                prev_expr                = call_sites.extra.bb.statements[prev_stmt_idx].expr
                cur_expr                 = ast
                live_in_before_prev_expr = call_sites.extra.live_in_before_prev_expr
                live_in_before_expr      = call_sites.extra.live_in_before_expr
                pattern_skeleton1        = (:(=), Tuple{SparseMatrixCSC,Int}, Expr(:call, TopNode(:indexed_next), Tuple{SparseMatrixCSC, SparseMatrixCSC}, Int, Int))
                pattern_skeleton2        = (:(=), SparseMatrixCSC, Expr(:call, TopNode(:getfield), :f1, 1))
                if match_skeletons(prev_expr, pattern_skeleton1,  call_sites, live_in_before_prev_expr, prev_expr, cur_expr, false) &&
                   match_skeletons(cur_expr, pattern_skeleton2,  call_sites, live_in_before_expr, prev_expr, cur_expr, false)
                    graph       = call_sites.extra.inter_dependence_graph
                    arr         = cur_expr.args[1]          # L
                    tup         = prev_expr.args[2].args[2] # GenSym(0)
                    index       = prev_expr.args[2].args[3] # x
                    tup_element = Expr(:TupleElement, tup, index)
                    build_vertices_and_edge(arr, tup_element, ROW_ROW, graph)
                    build_vertices_and_edge(arr, tup_element, COLUMN_COLUMN, graph)
                    return nothing
                end
            end
        
            rhs_type = type_of_ast_node(args[2], symbol_info)
            if !(rhs_type <: Vector) && !(rhs_type <: AbstractMatrix)
                # Look into rhs to see if its sub-expressions have anything to 
                # do with inter-dependence. For example, if rhs is dot(x, y),
                # although rhs type is not an array, x and y has ROW_ROW
                # inter-dependence 
                return build_inter_dependence_graph(args[2], call_sites)
            end
            lhs_type = type_of_ast_node(args[1], symbol_info)
            module_name   = ""
            function_name = ":="
            arg_types     = (lhs_type, rhs_type)
            arguments     = args
            return build_inter_dependence_graph_for_call(ast, module_name,
                                function_name, arg_types, arguments, call_sites)
        else
            throw(UnhandledExpr(head, args))
        end
    elseif asttyp == SymbolNode  || asttyp == Symbol   || asttyp == GenSym ||
           asttyp == LabelNode   || asttyp == GotoNode || asttyp == LineNumberNode ||
           asttyp == ASCIIString || asttyp == LambdaStaticData ||
           asttyp <: Number      || asttyp == NewvarNode || asttyp == GlobalRef ||
           asttyp == QuoteNode
        return nothing
    else
        throw(UnknownASTDistributivity(ast, typeof(ast)))
    end
    return nothing
end

@doc """ A handler to visit_expressions(). """
build_inter_dependence_graph(ast, call_sites :: CallSites, top_level_number, is_top_level, read) =
    build_inter_dependence_graph(ast, call_sites)

@doc """
Color the inter-dependence graph, starting with given vertex, which has been
colored itself.
"""
function color_inter_dependence_graph(
    liveness :: Liveness, 
    graph    :: InterDependenceGraph,
    from     :: InterDependenceGraphVertex,
    visited  :: Set{InterDependenceGraphVertex}
)
    if in(from, visited)
        return
    end    
    push!(visited, from)
    for (vertex, inverse) in from.neighbours
        color = inverse ? inverse_color_map[from.color] : from.color       
        if vertex.color != NOT_PERM_YET
            if vertex.color != color
                # Already colored, but in a different color. A conflict exists
                # But this can indicate that some colors are equal. For example,
                # from p = A * p, we know p.rowPermutation = A.rowPermutation=
                # A.colInversePermutation. Thus p's rowPermutation can have two
                # results (A.rowPermutation and A.colInversePermutation) but that
                # just means that A's rows and columns must be permuted in an
                # inverse relationship: P1 * A * P2, where P1 must equal P2'.
                # This is a nice constraint that compiler should tell libraries
                # so that rows and columns of A won't be reordered independently. 
                if vertex.color == NEVER_PERM
                    vertex2index = sort_inter_dependence_graph_vertices(graph)
                    dprintln(1, 0, "\nReordering cannot proceed: vertex ",
                             vertex2index[vertex], " in the inter-dependence ",
                             "graph cannot be reordered, as it might have array ",
                             "element access that the current implementation ", 
                             "does not support. However, its neighbour, vertex ",
                             vertex2index[from], " requires its color to be ",
                             permutation_color_to_str(color), ".")
                    dprintln(1, 0, "\nInter-dependence graph:")
                    dprintln(1, 1, liveness, graph)
                    throw(ConflictPermutation(from, vertex, color))
                else
                    dprintln(1, 0, "\nPermutation constraints found: ",
                             permutation_color_to_str(vertex.color),
                             " must equal ", 
                             permutation_color_to_str(color))
                    if (vertex.color == ROW_PERM && color == COL_INV_PERM) ||
                       (vertex.color == COL_INV_PERM && color == ROW_PERM) ||
                       (vertex.color == ROW_INV_PERM && color == COL_PERM) ||
                       (vertex.color == COL_PERM && color == ROW_INV_PERM)
                            graph.perm_restriction |= PERM_VECTORS_ARE_INVERSE
                    elseif (vertex.color == ROW_PERM && color == COL_PERM) ||
                           (vertex.color == COL_PERM && color == ROW_PERM) ||
                           (vertex.color == ROW_INV_PERM && color == COL_INV_PERM) ||
                           (vertex.color == COL_INV_PERM && color == ROW_INV_PERM)
                            graph.perm_restriction |= PERM_VECTORS_ARE_SAME
                    else
                        # There is no other known case so far
                        assert(false)
                    end
                end
            end
        else
            vertex.color = color
        end
        color_inter_dependence_graph(liveness, graph, vertex, visited)
    end
end

@doc """
Color the inter-dependence graph, starting with the seed in it.
"""
function color_inter_dependence_graph(
    liveness :: Liveness, 
    graph    :: InterDependenceGraph
)
    seed                      = graph.seed
    seed_rows_vertex          = graph.rows[seed]
    seed_columns_vertex       = graph.columns[seed]
    seed_rows_vertex.color    = ROW_PERM
    seed_columns_vertex.color = COL_PERM
    visited                   = Set{InterDependenceGraphVertex}()
    color_inter_dependence_graph(liveness, graph, seed_rows_vertex, visited)
    color_inter_dependence_graph(liveness, graph, seed_columns_vertex, visited)
end

@doc """
Add to the new expression a constant that indicates which permutation
vector to use for reordering
"""
function add_permutation_vector(
    new_expr :: Expr,
    perm     :: Int
)
    assert(perm == NOT_PERM_YET || perm == ROW_PERM || perm == ROW_INV_PERM || 
           perm == COL_PERM || perm == COL_INV_PERM)
    push!(new_expr.args, 
          perm == NOT_PERM_YET ? GlobalRef(Sparso, :NOT_PERM_YET) : 
          perm == ROW_PERM ? GlobalRef(Sparso, :ROW_PERM) :
          perm == ROW_INV_PERM ? GlobalRef(Sparso, :ROW_INV_PERM) :
          perm == COL_PERM ? GlobalRef(Sparso, :COL_PERM) :
          GlobalRef(Sparso, :COL_INV_PERM))
end

@doc """
Add to the new expression a matrix knob
"""
function add_matrix_knob(
    new_expr :: Expr,
    mknob    :: Symbol
)
    push!(new_expr.args, mknob)
end


@doc """
Add to the new expression an array for reordering or reverse reordering.
"""
function add_array(
    new_expr      :: Expr,
    A             :: Sym,
    symbol_info   :: Sym2TypeMap, 
    graph         :: InterDependenceGraph,
    matrices_done :: Bool,
    matrix2mknob  :: Dict{Symexpr, Symbol}
)
    if haskey(graph.rows, A) # If no key, the array does not appear in the loop at all
        if !matrices_done && type_of_ast_node(A, symbol_info) <: AbstractSparseMatrix
            Pr = graph.rows[A].color
            Pc = graph.columns[A].color
            if Pr != NOT_PERM_YET || Pc != NOT_PERM_YET
                push!(new_expr.args, A)
                add_permutation_vector(new_expr, Pr)
                add_permutation_vector(new_expr, Pc)
                add_matrix_knob(new_expr, matrix2mknob[A])
            end
        elseif  matrices_done && type_of_ast_node(A, symbol_info) <: Vector
            Pr = graph.rows[A].color
            if Pr != NOT_PERM_YET
                push!(new_expr.args, A)
                add_permutation_vector(new_expr, Pr)
            end
        end
    end
end

@doc """
Add to the new expression for reordering all the arrays live out of the 
statement that contains the reordering decider. 
This function should be called twice. The first time matrices_done = false, and 
the second time matrices_done = true. They will add matrices and vectors,
respectively.
"""
function add_arrays_to_reorder(
    new_expr      :: Expr,
    decider_stmt  :: Statement,
    symbol_info   :: Sym2TypeMap, 
    liveness      :: Liveness, 
    graph         :: InterDependenceGraph,
    FAR           :: Vector, # Vector{Symbol},
    matrices_done :: Bool,
    matrix2mknob  :: Dict{Symexpr, Symbol}
)
    live_out = LivenessAnalysis.live_out(decider_stmt, liveness)
    def      = LivenessAnalysis.def(decider_stmt, liveness)
#    use      = LivenessAnalysis.use(decider_stmt, liveness)
    
    # FAR (the symbols defined and used in the decider function call) have already
    # been reordered during the execution of that call. We expect a decide call
    # f happens in this way: 
    #   f(...)  or  x = f(...)
    # In the secon case, the defininition of x has also been reordered.
    # Reorder all the live-out symbols except them.
    # ISSUE: how to make sure f() is the only call happens in the statement it 
    # is called? What if for a statement like  x= f(...) * g(...)? 
    for A in setdiff(live_out, union(FAR, def))
        add_array(new_expr, A, symbol_info, graph, matrices_done, matrix2mknob)
    end
end

@doc """
Add to the new expression for reverse reodering all the arrays live out of 
from_bb and live into to_bb. 
This function should be called twice. The first time matrices_done = false, and 
the second time matrices_done = true. They will add matrices and vectors,
respectively.
"""
function add_arrays_to_reversely_reorder(
    new_expr      :: Expr,
    from_bb       :: BasicBlock,
    to_bb         :: BasicBlock,
    symbol_info   :: Sym2TypeMap, 
    liveness      :: Liveness, 
    graph         :: InterDependenceGraph,
    matrices_done :: Bool,
    matrix2mknob  :: Dict{Symexpr, Symbol}
)
    live_out = LivenessAnalysis.live_out(from_bb, liveness)
    live_in  = LivenessAnalysis.live_in(to_bb, liveness)
    for A in intersect(live_out, live_in)
        add_array(new_expr, A, symbol_info, graph, matrices_done, matrix2mknob)
    end
end

@doc """
Add to the new expression for reverse reordering all the arrays live into stmt.  
This function should be called twice. The first time matrices_done = false, and 
the second time matrices_done = true. They will add matrices and vectors,
respectively.
"""
function add_arrays_to_reversely_reorder(
    new_expr      :: Expr,
    stmt          :: Statement,
    symbol_info   :: Sym2TypeMap, 
    liveness      :: Liveness, 
    graph         :: InterDependenceGraph,
    matrices_done :: Bool,
    matrix2mknob  :: Dict{Symexpr, Symbol}
)
    live_in  = LivenessAnalysis.live_in(stmt, liveness)
    for A in live_in
        add_array(new_expr, A, symbol_info, graph, matrices_done, matrix2mknob)
    end
end

@doc """
Create reorder actions for the loop region.
"""
function create_reorder_actions(
    actions     :: Vector{Action},
    region      :: Region,
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    graph       :: InterDependenceGraph,
    cfg         :: CFG,
    call_sites  :: CallSites,
    FAR         :: Vector, # Vector{Symbol},
    fknob       :: Symbol,
)
    decider_bb       = call_sites.extra.decider_bb
    decider_stmt_idx = call_sites.extra.decider_stmt_idx
    matrix2mknob     = call_sites.extra.matrix2mknob
                               
    # Create an action that would insert new statements before the region
    # loop's head block (If it is a loop region), or at the beginning (if it is
    # a function region). There are two new statements: one is to set
    # the reordering decision maker, with the restrictioins of permutations,
    # and the other is to initialize reordering_status.
    region               = call_sites.region
    action_before_region = (typeof(region) == FunctionRegion) ? 
        InsertBeforeOrAfterStatement(Vector{Statement}(), true, region.entry, 1) :
        InsertBeforeLoopHead(Vector{Statement}(), region.loop, true)
    push!(actions, action_before_region)

    restriction = graph.perm_restriction
    stmt        = Expr(:call, GlobalRef(Sparso, :set_reordering_decision_maker), fknob, restriction)
    push!(action_before_region.new_stmts,  Statement(0, stmt))

    reordering_status = new_symbol("reordering_status")
    stmt              = Expr(:(=), reordering_status,
                              Expr(:call, TopNode(:vect), false, 
                                    GlobalRef(Main, :C_NULL), 
                                    GlobalRef(Main, :C_NULL),
                                    GlobalRef(Main, :C_NULL),
                                    GlobalRef(Main, :C_NULL),
                                    0.0))
    push!(action_before_region.new_stmts,  Statement(0, stmt))

    # In the region, after the statement that contains the reordering decision maker,
    # insert a statement to reorder other arrays. 
    inside_region_action = InsertBeforeOrAfterStatement(Vector{Statement}(),
                            false, decider_bb, decider_stmt_idx)
    push!(actions, inside_region_action)
    
    stmt = Expr(:call, GlobalRef(Sparso, :reordering),
                 fknob, reordering_status)
    push!(inside_region_action.new_stmts,  Statement(0, stmt))

    decider_stmt = decider_bb.statements[decider_stmt_idx]
    add_arrays_to_reorder(stmt, decider_stmt, symbol_info, liveness, graph, FAR, false, matrix2mknob)
    push!(stmt.args, QuoteNode(:__delimitor__))
    add_arrays_to_reorder(stmt, decider_stmt, symbol_info, liveness, graph, FAR, true, matrix2mknob)
    
    # At each region exit, insert reverse reordering of array.
    # Create statements that will delete the function knob at region exits
    if (typeof(region) == FunctionRegion)
        for pred_bb in region.exit.preds
            len = length(pred_bb.statements)
            assert(len > 0)
            last_stmt = pred_bb.statements[len]
            last_expr = last_stmt.expr
            assert(typeof(last_expr) <: Expr && last_expr.head == :return)

            action  = InsertBeforeOrAfterStatement(Vector{Statement}(), true, pred_bb, len)
            push!(actions, action)

            stmt = Expr(:call, GlobalRef(Sparso, :reverse_reordering),
                         reordering_status)
            push!(action.new_stmts,  Statement(0, stmt))
            
            add_arrays_to_reversely_reorder(stmt, last_stmt, symbol_info, liveness, graph, false, matrix2mknob)
            push!(stmt.args, QuoteNode(:__delimitor__))
            add_arrays_to_reversely_reorder(stmt, last_stmt, symbol_info, liveness, graph, true, matrix2mknob)
        end
    else
        for exit in region.exits
            action = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
            push!(actions, action)
    
            stmt = Expr(:call, GlobalRef(Sparso, :reverse_reordering),
                         reordering_status)
            push!(action.new_stmts,  Statement(0, stmt))
        
            add_arrays_to_reversely_reorder(stmt, exit.from_bb, exit.to_bb, symbol_info, liveness, graph, false, matrix2mknob)
            push!(stmt.args, QuoteNode(:__delimitor__))
            add_arrays_to_reversely_reorder(stmt, exit.from_bb, exit.to_bb, symbol_info, liveness, graph, true, matrix2mknob)
        end
    end
end

@doc """ 
Perform analyses for reordering. Write the intended transformation into actions.
"""
function reordering(
    actions     :: Vector{Action},
    region      :: Region,
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG,
    call_sites  :: CallSites
)
    old_extra   = call_sites.extra
    try
        fknob = call_sites.extra.reordering_decider_fknob
        if !in(fknob, call_sites.extra.function_knobs)
            # The AST that can decides reordering was probably removed due to,
            # for example, parent tree being changed.
            return actions
        end
        decider          = call_sites.extra.fknob2expr[fknob]
        FAR              = call_sites.extra.reordering_FAR
        seed             = FAR[1]
        matrix2mknob     = call_sites.extra.matrix2mknob
        call_sites.extra = ReorderingExtra(seed, decider, matrix2mknob)

        # An empty inter-dependence graph has been built. Build a row and a 
        # column vertex for the seed.
        graph = call_sites.extra.inter_dependence_graph
        build_vertex(seed, true, graph)
        build_vertex(seed, false, graph)
        
        recursive  = false
        visit_expressions(region, cfg, call_sites, recursive, build_inter_dependence_graph)

        color_inter_dependence_graph(liveness, graph)

        dprintln(1, 0, "\nColored inter-dependence graph:")
        dprintln(1, 1, liveness, graph)

        new_actions = Vector{Action}()
        create_reorder_actions(new_actions, region, symbol_info, liveness, graph, 
                               cfg, call_sites, FAR, fknob)

        dprintln(1, 0, "\nReordering actions to take:\n", new_actions)

        for action in new_actions
            push!(actions, action)
        end
    catch ex
        # In case any exception happen in the printing, try
        try
            Libc.flush_cstdio()
            flush(STDOUT)
            dprintln(1, 0, "Exception! Sparse Accelerator skips reordering the loop.")
            dprintln(1, 1, ex)
            dprintln(1, 0, "********************************************************************************")
            Libc.flush_cstdio()
            flush(STDOUT)
        catch
            # Do nothing
        end
    finally
        call_sites.extra = old_extra
        return actions
    end
end
