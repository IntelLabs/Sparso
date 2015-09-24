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
    
    InterDependenceGraphVertex(_array, _row_perm) = new(_array, _row_perm, NO_PERM,
                               Set{Tuple{InterDependenceGraphVertex, Bool}}())
end

@doc """
Inter-dependence graph. It separately represents rows and columns of arrays. It
will be colored starting from the seed.
"""
type InterDependenceGraph
    rows    :: Dict{Symexpr, InterDependenceGraphVertex}
    columns :: Dict{Symexpr, InterDependenceGraphVertex}
    seed    :: Sym
    
    InterDependenceGraph(_seed) = new(
        Dict{Symexpr, InterDependenceGraphVertex}(),
        Dict{Symexpr, InterDependenceGraphVertex}(),
        _seed)
end
 
@doc """
The CallSites' extra field for reordering.
"""
type ReorderingExtra
    seed                   :: Sym
    decider_ast            :: Expr
    decider_bb             :: Any # BasicBlock
    decider_stmt_idx       :: Int
    current_bb             :: Any # BasicBlock
    current_stmt_idx       :: Int
    inter_dependence_graph :: InterDependenceGraph

    ReorderingExtra(_seed, _decider_ast) = new(_seed, _decider_ast,
             nothing, 0, nothing, 0, InterDependenceGraph(_seed))
end

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
    array_index1 :: Int,
    array_index2 :: Int,
    relation     :: Int,
    ast          :: Expr,
    args         :: Vector,
    graph        :: InterDependenceGraph
)
    array1 = (array_index1 == 0) ? ast : args[array_index1]
    array2 = (array_index2 == 0) ? ast : args[array_index2]

    if relation == ROW_ROW
        vertex1 = build_vertex(array1, true, graph)
        vertex2 = build_vertex(array2, true, graph)
        build_edge(vertex1, vertex2, false)
    elseif relation == COLUMN_COLUMN
        vertex1 = build_vertex(array1, false, graph)
        vertex2 = build_vertex(array2, false, graph)
        build_edge(vertex1, vertex2, false)
    else
        assert(relation == COLUMN_ROW_INVERSE)
        vertex1 = build_vertex(array1, false, graph)
        vertex2 = build_vertex(array2, true, graph)
        build_edge(vertex1, vertex2, true)
    end
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
        call_sites.extra.decider_bb       = call_sites.extra.current_bb
        call_sites.extra.decider_stmt_idx = call_sites.extra.current_stmt_idx
    end
 
    if typeof(ast) == Expr
        symbol_info = call_sites.symbol_info
        head        = ast.head
        args        = ast.args

        module_name, function_name = "", ""
        if head == :call || head == :call1
            module_name, function_name = resolve_call_names(args)
            arg_types                  = ntuple(i-> type_of_ast_node(args[i+1],
                                                symbol_info), length(args) - 1)
            arguments                  = args[2 : end]
        elseif head == :(=)
            function_name = ":="
            arg_types     = ntuple(i-> type_of_ast_node(args[i],
                                   symbol_info), length(args))
            arguments     = args
        else
            throw(UnhandledExpr(head, args))
        end

        if function_name == ""
            throw(UnresolvedFunction(head, args[1]))
        end
        fd = look_for_function_description(module_name, function_name, arg_types)
        if fd == nothing
            throw(UndescribedFunction(module_name, function_name, arg_types))
        end
        if !fd.distributive
            throw(NonDistributiveFunction(module_name, function_name, arg_types))
        end
        for (array_index1, array_index2, relation) in fd.IA
            build_vertices_and_edge(array_index1, array_index2, relation, ast, arguments, graph)
        end
    end
    return nothing
end


@doc """ 
Build inter-dependence graph for the loop region.
"""
function build_inter_dependence_graph(
    region      :: LoopRegion,
    liveness    :: Liveness, 
    cfg         :: CFG,
    call_sites  :: CallSites
)
    blocks = cfg.basic_blocks
    for bb_index in region.loop.members
        bb                          = blocks[bb_index]
        statements                  = bb.statements
        call_sites.extra.current_bb = bb
        for stmt_idx in 1 : length(statements)
            call_sites.extra.current_stmt_idx = stmt_idx
            stmt                              = statements[stmt_idx]
            expr                              = stmt.expr
            if typeof(expr) != Expr
                continue
            end
            
            # Try to pattern match and replace this expression with ExprPatterns.
            CompilerTools.AstWalker.AstWalk(expr, build_inter_dependence_graph, call_sites)
        end
    end

    call_sites.extra.inter_dependence_graph
end

@doc """
Color the inter-dependence graph, starting with given vertex, which has been
colored itself.
"""
function color_inter_dependence_graph(
    graph :: InterDependenceGraph,
    from  :: InterDependenceGraphVertex
)
    assert(from.color != NO_PERM)
    for (vertex, inverse) in from.neighbours
        color = inverse ? inverse_color_map[from.color] : from.color
        if vertex.color != NO_PERM
            if vertex.color != color
                # Already colored, but in a different color. A conflict exists
                throw(ConflictPermutation(from, vertex, color))
            end
        else
            vertex.color = color
        end
        color_inter_dependence_graph(graph, vertex)
    end
end

@doc """
Color the inter-dependence graph, starting with the seed in it.
"""
function color_inter_dependence_graph(
    graph :: InterDependenceGraph
)
    seed                      = graph.seed
    seed_rows_vertex          = graph.rows[seed]
    seed_columns_vertex       = graph.columns[seed]
    seed_rows_vertex.color    = ROW_PERM
    seed_columns_vertex.color = COL_PERM
    color_inter_dependence_graph(graph, seed_rows_vertex)
    color_inter_dependence_graph(graph, seed_columns_vertex)
end

@doc """
Add to the new expression a constant that indicates which permutation
vector to use for reordering
"""
function add_permutation_vector(
    new_expr :: Expr,
    perm     :: Int
)
    assert(perm == NO_PERM || perm == ROW_PERM || perm == ROW_INV_PERM || 
           perm == COL_PERM || perm == COL_INV_PERM)
    push!(new_expr.args, 
          perm == NO_PERM ? GlobalRef(SparseAccelerator, :NO_PERM) : 
          perm == ROW_PERM ? GlobalRef(SparseAccelerator, :ROW_PERM) :
          perm == ROW_INV_PERM ? GlobalRef(SparseAccelerator, :ROW_INV_PERM) :
          perm == COL_PERM ? GlobalRef(SparseAccelerator, :COL_PERM) :
          GlobalRef(SparseAccelerator, :COL_INV_PERM))
end

@doc """
Add to the new expression an array for reordering or reverse reordering.
"""
function add_array(
    new_expr      :: Expr,
    A             :: Sym,
    symbol_info   :: Sym2TypeMap, 
    graph         :: InterDependenceGraph,
    matrices_done :: Bool
)
    if !matrices_done && type_of_ast_node(A, symbol_info) <: AbstractSparseMatrix
        Pr = graph.rows[A].color
        Pc = graph.columns[A].color
        if Pr != NO_PERM || Pc != NO_PERM
            push!(new_expr.args, A)
            add_permutation_vector(new_expr, Pr)
            add_permutation_vector(new_expr, Pc)
        end
    elseif  matrices_done && type_of_ast_node(A, symbol_info) <: Vector
        Pr = graph.rows[A].color
        if Pr != NO_PERM
            push!(new_expr.args, A)
            add_permutation_vector(new_expr, Pr)
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
    matrices_done :: Bool
)
    live_out = LivenessAnalysis.live_out(decider_stmt, liveness)
    for A in live_out
        add_array(new_expr, A, symbol_info, graph, matrices_done)
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
    liveness      :: Liveness, 
    graph         :: InterDependenceGraph,
    matrices_done :: Bool
)
    live_out = LivenessAnalysis.live_out(from_bb, liveness)
    live_in  = LivenessAnalysis.live_in(to_bb, liveness)
    for A in intersect(live_out, live_in)
        add_array(new_expr, A, symbol_info, graph, matrices_done)
    end
end

@doc """
Create reorder actions for the loop region.
"""
function create_reorder_actions(
    actions          :: Vector{Action},
    region           :: LoopRegion,
    symbol_info      :: Sym2TypeMap, 
    liveness         :: Liveness, 
    graph            :: InterDependenceGraph,
    cfg              :: CFG, 
    FAR              :: Vector{Symbol},
    fknob            :: Symbol,
    decider_bb       :: BasicBlock,
    decider_stmt_idx :: StatementIndex
)
    # Create an action that would insert new statements before the region
    # loop's head block. There are two new statements: one is to set
    # the reordering decision maker; the other is to initialize reordering_status.
    before_loop_action = InsertBeforeLoopHead(Vector{Statement}(), region.loop, true)
    push!(actions, before_loop_action)

    stmt = Expr(:call, GlobalRef(SparseAccelerator, :set_reordering_decision_maker), fknob)
    push!(before_loop_action.new_stmts,  Statement(0, stmt))

    reordering_status = gensym("reordering_status")
    stmt              = Expr(:(=), reordering_status,
                              Expr(:call, TopNode(:vect), false, 
                                    GlobalRef(Main, :C_NULL), 
                                    GlobalRef(Main, :C_NULL),
                                    GlobalRef(Main, :C_NULL),
                                    GlobalRef(Main, :C_NULL),
                                    0.0))
    push!(before_loop_action.new_stmts,  Statement(0, stmt))

    # In the loop, after the statement that contains the reordering decision maker,
    # insert a statement to reorder other arrays. 
    inside_loop_action = InsertBeforeOrAfterStatement(Vector{Statement}(),
                            false, decider_bb, decider_stmt_idx)
    push!(actions, inside_loop_action)
    
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :reordering),
                 fknob, reordering_status)
    push!(inside_loop_action.new_stmts,  Statement(0, stmt))

    add_arrays_to_reorder(stmt, liveness, false)
    push!(stmt.args, :__delimitor__)
    add_arrays_to_reorder(stmt, liveness, true)

    # At each loop exit, insert reverse reordering of array.
    # Create statements that will delete the function knob at region exits
    for exit in region.exits
        action = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
        push!(actions, action)

        stmt = Expr(:call, GlobalRef(SparseAccelerator, :reverse_reordering),
                       reordering_status)
        push!(action.new_stmts,  Statement(0, stmt))

        add_arrays_to_reversely_reorder(stmt, liveness, false)
        push!(stmt.args, :__delimitor__)
        add_arrays_to_reversely_reorder(stmt, liveness, true)
    end
end

@doc """ 
Perform analyses for reordering. Write the intended transformation into actions.
"""
function reordering(
    actions     :: Vector{Action},
    region      :: LoopRegion,
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG,
    call_sites  :: CallSites
)
    old_actions = copy(actions)
    old_extra   = call_sites.extra
    try
        decider = call_sites.extra.reordering_decider
        if decider == nothing
            return actions
        end

        FAR              = call_sites.extra.reordering_FAR
        seed             = FAR[1]
        fknob            = call_sites.extra.function_knobs[decider]
        call_sites.extra = ReorderingExtra(seed, decider)
        
        graph = build_inter_dependence_graph(region, liveness, cfg, call_sites)
        color_inter_dependence_graph(graph)

        dprintln(1, 0, "\nColored inter-dependence graph:")
        dprintln(1, 1, liveness, graph)

        create_reorder_actions(actions, region, symbol_info, liveness, graph, 
                               cfg, FAR, fknob, 
                               call_sites.extra.decider_bb,
                               call_sites.extra.decider_stmt_idx)

        dprintln(1, 0, "\nReordering actions to take:", actions)
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

        actions = old_actions
    finally
        call_sites.extra = old_extra
        return actions
    end
end