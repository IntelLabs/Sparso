# Reordering analysis: Feasibility, and where and what to reorder.

# ISSUES: 
# (1) We hard code a sparse matrix format to be SparseMatrixCSC{Cdouble, Cint},
# and a permutation and inverse permutation vector to be Vector{Cint}.
# Should make it general in future
# (2) We assume that there are NO aliases between any two arrays. We should
# use both static alias analysis and dynamic alias check to get rid of the 
# assumption.

@doc """ 
Create new statements that allocate space for permutation (and inverse 
permutation) vector. Return their symbols.
"""
function create_allocate_space_for_permutation(
    new_stmts :: Vector{Statement},
    M         :: Symbol
)
    P = gensym("P")
    stmt = :($P = Array(Cint, size($M, 2)))
    push!(new_stmts, Statement(0, stmt))

    inverse_P = gensym("inverse_P")
    stmt      = :($inverse_P = Array(Cint, size($M, 2)))
    push!(new_stmts, Statement(0, stmt))

    P, inverse_P
end

@doc """
Create new statements that will reorder matrix M with the permutation and inverse
permutation vector P and inverse_P, which will be computed if permute is true. 
"""
function create_reorder_matrix(
    new_stmts        :: Vector{Statement}, 
    M                :: Symbol, 
    P                :: Symbol, 
    inverse_P        :: Symbol, 
    permute          :: Bool,
    one_based_output :: Bool
)
    # Allocate space that stores the reordering result in Julia
    new_M = gensym(string(M))
    stmt  = Expr(:(=), new_M,
        Expr(:call, :SparseMatrixCSC, 
            Expr(:call, TopNode(:getfield), M, QuoteNode(:m)),
            Expr(:call, TopNode(:getfield), M, QuoteNode(:n)), 
            Expr(:call, :Array, :Cint, 
                Expr(:call, TopNode(:arraylen), 
                    Expr(:call, TopNode(:getfield), M, QuoteNode(:colptr)))),
            Expr(:call, :Array, :Cint, 
                Expr(:call, TopNode(:arraylen), 
                    Expr(:call, TopNode(:getfield), M, QuoteNode(:rowval)))),
            Expr(:call, :Array, :Cdouble, 
                Expr(:call, TopNode(:arraylen), 
                    Expr(:call, TopNode(:getfield), M, QuoteNode(:nzval))))
        )
    )
    push!(new_stmts, Statement(0, stmt))

    # Do the actual reordering through the C library
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :reorder_matrix),
                M, new_M, P, inverse_P, permute, one_based_output)
    push!(new_stmts, Statement(0, stmt))

    # Update the original matrix with the new data.
    stmt = Expr(:(=), M, new_M)
    push!(new_stmts, Statement(0, stmt))
end

@doc """
Create new statements that will reversely reorder matrix M with the permutation 
and inverse permutation vector P and inverse_P. 
"""
function create_reverse_reorder_matrix(
    new_stmts        :: Vector{Statement}, 
    M                :: Symbol, 
    P                :: Symbol, 
    inverse_P        :: Symbol, 
    one_based_output :: Bool
)
    create_reorder_matrix(new_stmts, M, inverse_P, P, false, one_based_output)
end

@doc """
Create new statements that will reorder vector V with the permutation vector P. 
"""
function create_reorder_vector(
    new_stmts :: Vector{Statement}, 
    V         :: Symbol, 
    P         :: Symbol
)
    # Allocate space that stores the reordering result in Julia
    new_V = gensym(string(V))
    stmt  = Expr(:(=), new_V,
                Expr(:call, :Array, :Cdouble, 
                    Expr(:call, TopNode(:arraylen), V)))
    push!(new_stmts, Statement(0, stmt))
    
    # Do the actual reordering through the C library
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :reorder_vector), V, new_V, P)
    push!(new_stmts, Statement(0, stmt))
    
    # Update the original vector with the new data
    stmt = Expr(:(=), V, new_V )
    push!(new_stmts, Statement(0, stmt))
end

@doc """
Create new statements that will reversely reorder vector V with the permutation 
vector P. 
"""
function create_reverse_reorder_vector(
    new_stmts :: Vector{Statement}, 
    V         :: Symbol, 
    P         :: Symbol
)
    # Allocate space that stores the reordering result in Julia
    new_V = gensym(string(V))
    stmt  = Expr(:(=), new_V,
                Expr(:call, :Array, :Cdouble, 
                    Expr(:call, TopNode(:arraylen), V)))
    push!(new_stmts, Statement(0, stmt))
    
    # Do the actual reordering through the C library
    stmt = Expr(:call, GlobalRef(SparseAccelerator, :reverse_reorder_vector), V, new_V, P)
    push!(new_stmts, Statement(0, stmt))
    
    # Update the original vector with the new data
    stmt = Expr(:(=), V, new_V )
    push!(new_stmts, Statement(0, stmt))
end

@doc """
Create an action that will insert new statements to reorder all arrays in X, 
with the permutation vector P and inverse permutation vector inverse_P.
"""
function create_reorder_X(
    new_stmts   :: Vector{Statement},
    X           :: Set,
    succ_bb     :: BasicBlock,
    P           :: Sym,
    inverse_P   :: Sym,
    symbol_info :: Sym2TypeMap
)
    for x in X
        if type_of_ast_node(x, symbol_info) <: AbstractMatrix
            create_reorder_matrix(new_stmts, x, P, inverse_P, true, true)
        elseif type_of_ast_node(x, symbol_info) <: AbstractVector
            create_reorder_vector(new_stmts, x, P)
        else
            throw(UnknownTypeToReorder(x, typeof(x, symbol_info)))
        end
    end
end

@doc """
Create an action that will insert new statements to reversely reorder all arrays 
in X, with the permutation vector P and inverse permutation vector inverse_P.
"""
function create_reverse_reorder_X(
    new_stmts   :: Vector{Statement},
    X           :: Set,
    succ_bb     :: BasicBlock,
    P           :: Sym,
    inverse_P   :: Sym,
    symbol_info :: Sym2TypeMap
)
    for x in X
        if type_of_ast_node(x, symbol_info) <: AbstractMatrix
            create_reverse_reorder_matrix(new_stmts, x, P, inverse_P, true)
        elseif type_of_ast_node(x, symbol_info) <: AbstractVector
            create_reverse_reorder_vector(new_stmts, x, P)
        else
            throw(UnknownTypeToReorder(x, typeof(x, symbol_info)))
        end
    end
end

@doc """
Find a SpMV from a call site.
"""
function discover_a_SpMV(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    if typeof(ast) <: Expr
        head = ast.head
        if head == :call || head == :call1
            args = ast.args
            module_name, function_name = resolve_call_names(args)
            if module_name == "Main" && function_name == "*" &&
                length(args) == 3 && 
                type_of_ast_node(args[2], call_sites.symbol_info) <: SparseMatrixCSC &&
                type_of_ast_node(args[3], call_sites.symbol_info) <: Vector
                    site = CallSite(ast)
                    push!(call_sites.sites, site) 
            end
        end
    end
    return nothing
end

@doc """
Find all SpMVs from the loop region.
"""
function discover_SpMVs(
    actions       :: Vector{Action},
    region        :: LoopRegion,
    symbol_info   :: Sym2TypeMap,
    cfg           :: CFG
)
    L      = region.loop
    blocks = cfg.basic_blocks
    SpMVs  = CallSites(Set{CallSite}(), region, symbol_info, 
                       Vector{Pattern}(), actions, nothing)
    for bb_idx in L.members
        bb         = blocks[bb_idx]
        statements = bb.statements
        for stmt_idx in 1 : length(statements)
            stmt = statements[stmt_index]
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            CompilerTools.AstWalker.AstWalk(expr, discover_a_SpMV, SpMVs)
        end
    end
    SpMVs
end

@doc """
For some SpMV calls, replace them with SpMV calls that also measure the benefit
of reordering. So far, we do this only for the first SpMV call.
TODO: find out the most time-consuming SpMVs, and decides the benefit from them. 
"""
function SpMVs_with_reordering_benefit(
    SpMVs              :: CallSites,
    first_reorder_done :: Symbol, 
    beneficial         :: Symbol
)
    for site in SpMVs.sites
        func = LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), 
            :SparseAccelerator, QuoteNode(:SpMV_conditional_reordering))
        site.ast.args = [func; site.ast.args[2:3]; first_reorder_done; beneficial]
        break
    end
end

@doc """
Create reorder actions for the loop region based on the IN/OUT of the nodes in
the reorder graph.
"""
function create_reorder_actions(
    actions       :: Vector{Action},
    reorder_graph :: ReorderGraph,
    region        :: LoopRegion,
    symbol_info   :: Sym2TypeMap, 
    liveness      :: Liveness, 
    cfg           :: CFG, 
    FAR           :: Vector{Symbol}
)
    loop_head = region.loop.head
    
    # Test if we should do reordering conditionally, i.e. based on a cost-
    # benefit model. So far, the model is only based on the existence of a
    # SpMV inside the loop; for that SpMV, we will let it run once, measure
    # how far its memory bandwidth is from the machine's peak memory bandwidth;
    # if it is far, reordering is potentially good, and we turn it on.
    beneficial             = gensym("beneficial")
    first_reorder_done     = gensym("first_reorder_done")
    conditional_reordering = false
    if reorder_when_beneficial
        SpMVs = discover_SpMVs(actions, lives, loop_info, symbol_info, region)
        if !isempty(SpMVs.sites)
            # TODO: turn on conditional reordering only when SpMVs are dominating
            # the execution time of the loop region, for example, when there is
            # no other more expensive matrix operations than SpMVs
            conditional_reordering = true
            SpMVs_with_reordering_benefit(SpMVs, first_reorder_done, beneficial)
        end            
    end

    if conditional_reordering
        # Insert new statements before the loop to initialize the conditions
        init_action = InsertBeforeLoopHead(Vector{Statement}(), region.loop, true)
        push!(actions, init_action)
        
        init_stmt = Statement(0, 
            Expr(:call, GlobalRef(SparseAccelerator, :init_conditional_reordering), 
                first_reorder_done, beneficial))
        push!(init_action.new_stmts, init_stmt)
    end

    # Find the first sparse matrix in FAR. Use it to get permutation array
    M = nothing
    for M in FAR
        if type_of_ast_node(M, symbol_info) <: AbstractSparseMatrix
            break
        end
    end
    assert(M !=nothing)

    # Create an action that would insert new statements before the region
    # loop's head block. In case of conditional reordering, it is before
    # the head, but inside the loop
    outside_loop         = !conditional_reordering
    first_reorder_action = InsertBeforeLoopHead(Vector{Statement}(), region.loop, outside_loop)
    push!(actions, first_reorder_action)
    
    if conditional_reordering
        # Check if reordering is not turned on, and it is beneficial to do so
        goto_stmt = Statement(0, Expr(:gotoifnot, 
            Expr(:&&, Expr(:call, TopNode(:arrayref), beneficial, 1), 
                Expr(:call, :!, first_reorder_done)), loop_head.label))
        push!(first_reorder_action.new_stmts, goto_stmt)
        
        done_stmt = Statement(0, Expr(:(=), first_reorder_done, true))
        push!(first_reorder_action.new_stmts, done_stmt)
    end

    # Allocate space to store the permutation and inverse permutation info
    P, inverse_P = create_allocate_space_for_permutation(first_reorder_action.new_stmts, M)

    # Compute P and inverse_P, and reorder M
    create_reorder_matrix(first_reorder_action.new_stmts, M, P, inverse_P, true, true)

    # For every edge pred->succ, for any x live into succ, if 
    # (1) x is in pred.Out but not in succ.In, create action "reverse_reorder(x)"
    # (2) x is not in pred.Out but in succ.In, create action "reorder(x)"
    # Insert the action on the edge between pred and succ.
    blocks = cfg.basic_blocks
    vertices = union(reorder_graph.vertices_in_region, reorder_graph.vertices_outside_region)
    for succ in vertices
        if succ.bb_idx == PSEUDO_BLOCK_INDEX
            continue
        end

        succ_bb         = blocks[succ.bb_idx]
        live_in_succ_bb = LivenessAnalysis.live_in(succ_bb, liveness)
        for pred in succ.preds
            if pred.kind == RG_NODE_ENTRY
                # Special handling for entry.
                # In the current implementation, we have only one entry for only
                # one kind of region (loop region), and it has one and only one 
                # successor: the loop head. Avoid reordering M again, which has
                # been reordered above.
                assert(length(pred.succs) == 1)
                assert(succ.bb_idx == loop_head)
                X = intersect(pred.Out, live_in_succ_bb)
                m = Set(); push!(m, M)
                setdiff!(X, m)
                create_reorder_X(first_reorder_action.new_stmts, X, succ_bb, P, inverse_P, symbol_info)
            else
                pred_bb = blocks[pred.bb_idx]
                X       = setdiff(succ.In, pred.Out)
                X       = intersect(X, live_in_succ_bb)
                X1      = setdiff(pred.Out, succ.In)
                X1      = intersect(X1, live_in_succ_bb)
                if !isempty(X) || !isempty(X1)
                    action  = InsertOnEdge(Vector{Statement}(), pred_bb, succ_bb)
                    push!(actions, action)
                    if conditional_reordering
                        # Check if reordering has been turned on
                        check_stmt = Expr(:gotoifnot, first_reorder_done, succ.bb_idx)
                        push!(actions.new_stmts, check_stmt)
                    end
                    if !isempty(X)
                        create_reorder_X(action.new_stmts, X, succ_bb, P, inverse_P, symbol_info)
                    end
                    if !isempty(X1)
                        create_reverse_reorder_X(action.new_stmts, X1, succ_bb, P, inverse_P, symbol_info)
                    end
                end
            end
        end
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
    if call_sites.extra.reordering_decider == nothing
        return actions
    end
    assert(!isempty(call_sites.extra.reordering_FAR))

    decider = call_sites.extra.reordering_decider
    FAR     = call_sites.extra.reordering_FAR
    fknob   = call_sites.extra.function_knobs[decider]

    distributive = check_distributivity(region, cfg, symbol_info)
    if distributive
        stmt_clusters = find_inter_dependent_arrays(region, cfg, symbol_info)
        reorder_graph = discover_reorderable_arrays(region, stmt_clusters, liveness, cfg, decider, FAR, fknob)
        create_reorder_actions(actions, reorder_graph, region, symbol_info, liveness, cfg, FAR)
    end

    dprintln(1, 0, "\nReordering actions to take:", actions)
    
    actions
end

# The functions below are obsolete, as we change to let the library to decide
# reordering or not, permutation vectors, etc.

@doc """ 
Find the first arrays to reorder from the function's input parameters. So far,
 the first sparse matrix paramter is regarded as the only first array to reorder.
"""
function find_first_arrays_to_reorder(
    func_ast    :: Expr, 
    symbol_info :: Sym2TypeMap
)
    assert(func_ast.head == :lambda)
    lambda = lambdaExprToLambdaInfo(func_ast)
    FAR    = Vector{Symbol}()
    for i in lambda.input_params
        if type_of_ast_node(i, symbol_info) <: AbstractSparseMatrix
            push!(FAR, i)
            break
        end
    end
    return FAR
end

@doc """ 
Perform analyses for reordering. Write the intended transformation into actions.
"""
function reordering(
    actions     :: Vector{Action},
    regions     :: Vector{LoopRegion},
    func_ast    :: Expr, 
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG, 
    loop_info   :: DomLoops
)
    FAR     = find_first_arrays_to_reorder(func_ast, symbol_info)
    if isempty(FAR)
        return actions
    end
    
    for region in regions
        distributive = check_distributivity(region, cfg, symbol_info)
        if distributive
            stmt_clusters = find_inter_dependent_arrays(region, cfg, symbol_info)
            reorder_graph = discover_reorderable_arrays(region, stmt_clusters, liveness, cfg, FAR)
            create_reorder_actions(actions, reorder_graph, region, symbol_info, liveness, cfg, FAR)
        end
    end

    dprintln(1, 0, "\nReordering actions to take:", actions)
    
    actions
end