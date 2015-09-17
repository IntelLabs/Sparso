# Reorderable array discovery analysis

@doc """
A vertex in a reorder graph, which is for reorderable array discovery analysis 
for a region. 

There are vertices inside and outside the region. 
(1) Inside the region: 
    A NORMAL vertex corresponds to one statement in a basic block in the region. 
    An EMPTY vertex is inserted for an empty block in the region. 
    An ENTRY is a special vertex inserted for the region, representing an entry
    boundary of the region.
(2) Outside the region:
    An OUTSIDE vertex is created for a predecessor or successor of a vertex, 
    when the predecessor/successor is outside the region. 

The difference between vertices inside and outside the region is: those inside
the region are subject to dataflow propagation, while those outside are not.
The outside vertices only provides a fixed IN/OUT (empty set) to its 
predecessor/successor vertices inside.
"""
const RG_NODE_ENTRY   = 0xfff0
const RG_NODE_EMPTY   = 0xfff1
const RG_NODE_NORMAL  = 0xfff2
const RG_NODE_OUTSIDE = 0xfff4

const PSEUDO_BLOCK_INDEX     = -0xfffffff
const PSEUDO_STATEMENT_INDEX = -0xfffffff

type ReorderGraphVertex
    bb_idx   :: BasicBlockIndex
    stmt_idx :: StatementIndex
    kind     :: Int # One of RG_NODE_ENTRY/EMPTY/NORMAL/OUTSIDE
    succs    :: Set{ReorderGraphVertex}
    preds    :: Set{ReorderGraphVertex}
    In       :: Set # Set{Sym} causes incompatibility with Liveness package's TopLevelStatement.def and use
    Out      :: Set

    ReorderGraphVertex(basic_block_index, statement_idx, node_kind) = 
        new(basic_block_index, statement_idx, node_kind,
            Set{ReorderGraphVertex}(), Set{ReorderGraphVertex}(), Set(), Set())
end

@doc """
Reorder graph for a region, used for discovering reorderable arrays. It contains
vertices in and outside the region.
"""
type ReorderGraph
    vertices_in_region      :: Set{ReorderGraphVertex}
    vertices_outside_region :: Set{ReorderGraphVertex}
    stmt_clusters           :: Statement2Clusters      # Only for debugging purpose
    
    ReorderGraph(_stmt_clusters) = new(Set{ReorderGraphVertex}(), Set{ReorderGraphVertex}(), _stmt_clusters)
end

function build_reorder_graph_for_region(
    region        :: LoopRegion,
    stmt_clusters :: Statement2Clusters,
    liveness      :: Liveness, 
    cfg           :: CFG
)
    graph = ReorderGraph(stmt_clusters)
    
    L = region.loop
    if isempty(L.members) 
        return graph
    end

    # Build vertices and connect vertices in the same block
    first_node = Dict{BasicBlockIndex, ReorderGraphVertex}()
    last_node  = Dict{BasicBlockIndex, ReorderGraphVertex}()
    blocks     = cfg.basic_blocks
    for bb_idx in L.members
        bb = blocks[bb_idx]
        
        # Create vertices (inside the region) for the statements in the block
        # And connect them
        statements = bb.statements
        if length(statements) == 0
            # Empty block. Make a pseudo vertex
            vertex = ReorderGraphVertex(bb_idx, stmt_idx, RG_NODE_EMPTY)
            push!(graph.vertices_in_region, vertex)
            first_node[bb_idx] = last_node[bb_idx] = vertex
        else 
            prev = nothing
            for stmt_idx in 1 : length(statements)
                vertex = ReorderGraphVertex(bb_idx, stmt_idx, RG_NODE_NORMAL)
                push!(graph.vertices_in_region, vertex)
                
                if stmt_idx == 1
                    first_node[bb_idx] = vertex
                end
                if stmt_idx == length(statements)
                    last_node[bb_idx] = vertex
                end
                
                if prev != nothing
                    push!(prev.succs, vertex)
                    push!(vertex.preds, prev)
                end
                prev = vertex
            end
        end

        # Create vertices outside the region if predecessor/successor blocks are
        # outside the region.
        for pred in bb.preds
            pred_bb_idx = pred.label
            if !in(pred_bb_idx, L.members) && !haskey(first_node, pred_bb_idx)
                # if pred is the predecessor of the loop head outside the loop,
                # Build a pseudo entry vertex
                if bb_idx == L.head
                    vertex = ReorderGraphVertex(pred_bb_idx, PSEUDO_STATEMENT_INDEX, RG_NODE_ENTRY)
                    push!(graph.vertices_in_region, vertex)
                else
                    vertex = ReorderGraphVertex(pred_bb_idx, PSEUDO_STATEMENT_INDEX, RG_NODE_OUTSIDE)
                    push!(graph.vertices_outside_region, vertex)
                end
                first_node[pred_bb_idx] = last_node[pred_bb_idx] = vertex
            end
        end
        for succ in bb.succs
            succ_bb_idx = succ.label
            if !in(succ_bb_idx, L.members) && !haskey(first_node, succ_bb_idx)
                vertex = ReorderGraphVertex(succ_bb_idx, PSEUDO_STATEMENT_INDEX, RG_NODE_OUTSIDE)
                push!(graph.vertices_outside_region, vertex)
                first_node[succ_bb_idx] = last_node[succ_bb_idx] = vertex
            end
        end
    end

    # Connect first/last vertices from different blocks
    for bb_idx in L.members
        bb = blocks[bb_idx]
        
        for pred in bb.preds
            pred_bb_idx = pred.label
            push!(first_node[bb_idx].preds, last_node[pred_bb_idx])
            push!(last_node[pred_bb_idx].succs, first_node[bb_idx])
        end

        for succ in bb.succs
            succ_bb_idx = succ.label
            push!(first_node[succ_bb_idx].preds, last_node[bb_idx])
            push!(last_node[bb_idx].succs, first_node[succ_bb_idx])
        end
    end

    return graph
end

const UNIVERSE_SYM = gensym("universe")

@doc """
Intersection of all the predecessor vertices' OUT.
"""
function intersect_preds_out(vertex :: ReorderGraphVertex)
    X = Set()
    first = true
    for pred in vertex.preds
        if first || in(UNIVERSE_SYM, X)
            X = copy(pred.Out)
            first = false
        else
            if !in(UNIVERSE_SYM, pred.Out)
                X = intersect(X, pred.Out)
            end 
        end
    end
    X
end

@doc """
Intersection of all the successor vertices' IN.
"""
function intersect_succs_in(vertex)
    X = Set()
    first = true
    for succ in vertex.succs
        if first || in(UNIVERSE_SYM, X)
            X = copy(succ.In)
            first = false
        else
            if !in(UNIVERSE_SYM, succ.In)
                X = intersect(X, succ.In)
            end 
        end
    end
    X
end

@doc """
Forward transfer of input X through vertex B.
"""
function forward_transfer(
    B             :: ReorderGraphVertex, 
    X             :: Set, 
    stmt_clusters :: Statement2Clusters, 
    cfg           :: CFG
)
    assert(B.kind != RG_NODE_OUTSIDE)
    if B.kind == RG_NODE_ENTRY || B.kind == RG_NODE_EMPTY
        return copy(X)
    end
    
    bb       = cfg.basic_blocks[B.bb_idx]
    stmt     = bb.statements[B.stmt_idx]
    clusters = stmt_clusters[stmt]
    if isempty(clusters)
        return copy(X)
    end
    
    result = Set()
    for x in X
        for cluster in clusters
            if in(x, cluster.RHS)
                union!(result, cluster.RHS)
                # no matter a RHS symbol is in LHS or not, it should be transfered:
                # If it is in, then of course it will (It means the new def of the symbol)
                # If not, it will as well (It means the current def of the symbol)
                union!(result, cluster.LHS)
            else
                if in(x, cluster.LHS)
                    # The current def of the symbol is killed. The new def, unless
                    # some symbol in the RHS is inter-dependent on it, will not be added
                    # Do nothing.
                else
                    push!(result, x)
                end
            end
        end
    end
    return result
end

@doc """
Backward transfer of input X through vertex B.
"""
function backward_transfer(
    B             :: ReorderGraphVertex, 
    X             :: Set, 
    stmt_clusters :: Statement2Clusters, 
    cfg           :: CFG
)
    assert(B.kind != RG_NODE_OUTSIDE)
    if B.kind == RG_NODE_ENTRY || B.kind == RG_NODE_EMPTY
        return copy(X)
    end
    
    bb       = cfg.basic_blocks[B.bb_idx]
    stmt     = bb.statements[B.stmt_idx]
    clusters = stmt_clusters[stmt]
    if isempty(clusters)
        return copy(X)
    end
    
    result = Set()
    for x in X
        for cluster in clusters
            if in(x, cluster.LHS)
                union!(result, cluster.RHS)
            else
                if in(x, cluster.RHS)
                    union!(result, cluster.RHS)
                else
                    push!(result, x)
                end
            end
        end
    end
    return result
end

@doc """
A forward pass of dataflow propagation. Preconditioning is specially handled.
Return true if IN or OUT of any vertex has ever changed during the propagation.
"""
function forward_pass(
    vertices        :: Set{ReorderGraphVertex}, 
    stmt_clusters   :: Statement2Clusters, 
    cfg             :: CFG,
    preconditioning :: Bool
)
    ever_changed = false
    changed      = true
    while changed
        changed = false
        for vertex in vertices
            if preconditioning && vertex.kind == RG_NODE_ENTRY
                continue
            end
            X = intersect_preds_out(vertex)
            if !preconditioning
                X = union(vertex.In, X)
            end
            if !(X == vertex.In)
                changed      = true
                ever_changed = true
                vertex.In    = X
            end
            X = forward_transfer(vertex, vertex.In, stmt_clusters, cfg)
            if !preconditioning
                X = union(vertex.Out, X)
            end
            if !(X == vertex.Out)
                changed      = true
                ever_changed = true
                vertex.Out   = X
            end
        end
    end
    ever_changed
end

@doc """
A backward pass of dataflow propagation.
Return true if IN or OUT of any vertex has ever changed during the propagation.
"""
function backward_pass(
    vertices        :: Set{ReorderGraphVertex}, 
    stmt_clusters   :: Statement2Clusters, 
    cfg             :: CFG
)
    ever_changed = false
    changed      = true
    while changed
        changed = false
        for vertex in vertices
            X = intersect_succs_in(vertex)
            X = union(vertex.Out, X)
            if !(X == vertex.Out)
                changed      = true
                ever_changed = true
                vertex.Out   = X
            end
            X = backward_transfer(vertex, vertex.Out, stmt_clusters, cfg)
            X = union(vertex.In, X)
            if !(X == vertex.In)
                changed      = true
                ever_changed = true
                vertex.In    = X
            end
        end
    end
    ever_changed
end

@doc """
Find what arrays to be reordered and where.
"""
function discover_reorderable_arrays(
    region          :: Region, 
    stmt_clusters   :: Statement2Clusters, 
    liveness        :: Liveness, 
    cfg             :: CFG,
    FAR             :: Vector{Symbol}
)
    reorder_graph = build_reorder_graph_for_region(region, stmt_clusters, liveness, cfg)

    dprintln(1, 0, "Initial reorder graph:")
    dprintln(1, 1, liveness, reorder_graph)

    # Do bi-directional dataflow analysis on the graph.
    vertices_in_region = reorder_graph.vertices_in_region
    
    # Step 1: initialization
    for vertex in vertices_in_region
        if vertex.kind == RG_NODE_ENTRY
            union!(vertex.Out, FAR)
         else
            push!(vertex.Out, UNIVERSE_SYM)
        end
    end
    dprintln(1, 0, "Initialized reorder graph:")
    dprintln(1, 1, liveness, reorder_graph)

    # Step 2: preconditioning
    forward_pass(vertices_in_region, stmt_clusters, cfg, true)

    dprintln(1, 0, "Preconditioned reorder graph:")
    dprintln(1, 1, liveness, reorder_graph)

    # Step 3: Growth (repetitive backward and forward pass)
    changed = true
    i       = 1
    while changed
        changed  = backward_pass(vertices_in_region, stmt_clusters, cfg)
        
        dprintln(1, 0, "Reorder graph after backward pass ", i)
        dprintln(1, 1, liveness, reorder_graph)

        changed |= forward_pass(vertices_in_region, stmt_clusters, cfg, false)

        dprintln(1, 0, "Reorder graph after forward pass ", i)
        dprintln(1, 1, liveness, reorder_graph)
        
        i = i + 1
    end
        
    return reorder_graph
end