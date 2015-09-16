# This file is for inter-dependent array analysis.

@doc """ 
A set of inter-dependent arrays in the same statement. Some of them appear 
in the Left Hand Side (LHS) of the statement, some in the Right Hand side(RHS).
"""
type Cluster
    LHS :: Set{Sym}
    RHS :: Set{Sym}
    
    Cluster() = new(Set{Sym}(), Set{Sym}())
end

@doc """ 
A map from a statement to all the clusters of it. Each cluster contains 
some arrays that are inter-dependent with each other.
"""
typealias Statement2Clusters Dict{Statement, Set{Cluster}}

@doc """ 
The expression tree starting with a root, located in either left or right hand
side of a statement. 
"""
type ExpressionTree
    root :: Any
    left :: Bool
end

@doc """ 
If the node is an array, and is a terminal in its AST, add it to the cluster.
"""
function add_terminal_array_to_cluster(
    node        :: Any, 
    left        :: Bool, 
    cluster     :: Cluster, 
    symbol_info :: Sym2TypeMap
)
    if type_of_ast_node(node, symbol_info) <: AbstractArray
        # Symbol, SymbolNode and GenSym are terminals of the AST.
        if typeof(node) == Symbol
            push!(left ? cluster.LHS : cluster.RHS, node)
        elseif typeof(node) == SymbolNode
            push!(left ? cluster.LHS : cluster.RHS, node.name)
        elseif typeof(node) == GenSym
            push!(left ? cluster.LHS : cluster.RHS, node.id)
        end
    end 
end

@doc """ 
Starting from the AST node, which is either in the left or right hand side of
a statement, add all the arrays at the terminals of the AST into the cluster. 
Break the AST into sub-trees when a scalar is encountered.
"""
function find_inter_dependent_arrays(
    node        :: Any, 
    left        :: Bool,
    cluster     :: Cluster,
    trees       :: Set{ExpressionTree},
    symbol_info :: Sym2TypeMap
)
    add_terminal_array_to_cluster(node, left, cluster, symbol_info)
    
    if typeof(node) <: Expr
        head, args = node.head, node.args
        if head == :call || head == :call1
            # Compose a type tuple for the arguments, except he first argument, 
            # which is the function.
            arg_types = ntuple(i-> type_of_ast_node(args[i+1], symbol_info), length(args) - 1)
            all_numbers, some_arrays = numbers_or_arrays(node.typ, arg_types)
            if all_numbers || !some_arrays
                # The function call's result and arguments are all numbers, or 
                # some are non-numbers (like Range{UInt64}) but not regular arrays. 
                # No arrays to care for. Do nothing
            else
                module_name, function_name = resolve_call_names(args)
                if function_name == ""
                    throw(UnresolvedFunction(head, args[1]))
                end
                fd = look_for_function_description(module_name, function_name, arg_types)
                if fd != nothing
                    for S in fd.IA
                        for x in S
                            if x != 0 
                                # Not the result of the function call. Continue
                                # to further decompose/cluster from it
                                find_inter_dependent_arrays(node.args[x + 1], left, cluster, trees, symbol_info)
                            end
                        end
                    end
                else
                    throw(UndescribedFunction(module_name, function_name, arg_types))
                end
            end
        end
    
        # Recusively visit each argument to further decompose/cluster
        for c in node.args
            if type_of_ast_node(c, symbol_info) <: Number
                # Decompose: break the tree from this node as another tree. 
                # The arrays inside that tree have nothing to do with those in 
                # this tree.
                push!(trees, ExpressionTree(c, left))
            else
                find_inter_dependent_arrays(c, left, cluster, trees, symbol_info)
            end
        end
    end
end

@doc """ 
For each statement in the region, decompose it into some expression trees.
The arrays at the terminals of the same expression tree are inter-dependent 
with each other, and thus compose of a cluster.
"""
function find_inter_dependent_arrays(
    region      :: LoopRegion,
    cfg         :: CFG, 
    symbol_info :: Sym2TypeMap
)
    mapping = Statement2Clusters()
    blocks  = cfg.basic_blocks
    for bb_index in region.loop.members
        bb = blocks[bb_index]
        for stmt in bb.statements
            mapping[stmt] = Set{Cluster}()
            trees         = Set{ExpressionTree}()
            push!(trees, ExpressionTree(stmt.expr, false))
            while !isempty(trees)
                tree    = pop!(trees)
                cluster = Cluster()
                if typeof(tree.root) == Expr && tree.root.head == :(=)
                    # An assignment can happen only at top level, not inside LHS or RHS 
                    assert(tree.root == stmt.expr)
                    find_inter_dependent_arrays(tree.root.args[1], true, cluster, trees, symbol_info)
                    find_inter_dependent_arrays(tree.root.args[2], false, cluster, trees, symbol_info)
                else
                    find_inter_dependent_arrays(tree.root, tree.left, cluster, trees, symbol_info)
                end
                if !isempty(cluster.LHS) || !isempty(cluster.RHS)
                    push!(mapping[stmt], cluster)
                end
            end
            
        end

    end
    return mapping
end
