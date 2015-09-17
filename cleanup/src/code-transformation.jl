@doc """
Group InsertBeforeOrAfterStatement actions in decreasing order of
(basic block index, statement index) of the actions; for the same
statement, "after" is before "before".
"""
function group_action(
    insert_before_or_after_statement_actions :: Vector{InsertBeforeOrAfterStatement},
    action                                   :: InsertBeforeOrAfterStatement
)
    for i in 1 : length(insert_before_or_after_statement_actions)
        action1 = insert_before_or_after_statement_actions[i]
        if action.bb == action1.bb 
            if action.stmt_idx == action1.stmt_idx
                if !action.before && action1.before
                    # Let the actions for a later statement happen earlier since
                    # they do not impact the positions of earlier statements.
                    insert!(insert_before_or_after_statement_actions, i, action)
                    
                    # No need to try match with other elements: the vector is
                    # ensured to have only one elment for one "before" and/or
                    # one "after" for one specific statement.
                    return
                elseif action.before && !action1.before
                    continue
                else
                    # The two actions are both before or both after the statement.
                    append!(action1.new_stmts, action.new_stmts)
                    return
                end
            elseif action.stmt_idx > action1.stmt_idx
                # Let the actions for a later statement happen earlier since
                # they do not impact the positions of earlier statements.
                insert!(insert_before_or_after_statement_actions, i, action)
                return
            end
        elseif action.bb.label > action.bb.label
            insert!(insert_before_or_after_statement_actions, i, action)
            return
        end
    end
    push!(insert_before_or_after_statement_actions, action)
end

@doc """
Group InsertBeforeLoopHead actions.
"""
function group_action(
    insert_before_loop_head_actions :: Vector{InsertBeforeLoopHead},
    action                          :: InsertBeforeLoopHead
)
    for action1 in insert_before_loop_head_actions
        if action1.loop == action.loop && action1.outside_loop == action.outside_loop
            append!(action1.new_stmts, action.new_stmts)
            
            # No need to try match with other elements: the vector is ensured
            # to have only one elment for one specific block.
            return
        end
    end
    push!(insert_before_loop_head_actions, action)
end

@doc """
Group InsertOnEdge actions.
"""
function group_action(
    insert_on_edge_actions :: Vector{InsertOnEdge},
    action                 :: InsertOnEdge
)
    for action1 in insert_on_edge_actions
        if action1.from_bb == action.from_bb && action1.to_bb == action.to_bb
            append!(action1.new_stmts, action.new_stmts)
            
            # No need to try match with other elements: the vector is ensured
            # to have only one elment for one specific edge.
            return
        end
    end
    push!(insert_on_edge_actions, action)
end

@doc """
The index of the statement in the block. Return -1 if the statement is not in 
the block.
"""
function index_of_statement(
    bb   :: BasicBlock,
    stmt :: Statement
)
    for idx = 1 : length(bb.statements)
        if bb.statements[idx] == stmt
            return idx
        end
    end
    return -1
end

@doc """
Perform all the actions. Return the transformed AST.
"""
function code_transformation(
    actions  :: Vector{Action},
    cfg      :: CFG
)
    # Implement each action. Some important details to care:
    # (1) InsertBeforeLoopHead actions must be done last, because it changes the
    # predecessor-successor relationship, which affects InsertOnEdge (Some edge
    # no longer exists after InsertBeforeLoopHead).
    # (2) Actions that work on the same loop for InsertBeforeLoopHead, or the same
    # Statement for InsertBeforeOrAfterStatement, or the same edge for InsertOnEdge, 
    # should be grouped and implemented together
    
    # Group the actions.
    insert_before_or_after_statement_actions = Vector{InsertBeforeOrAfterStatement}()
    insert_before_loop_head_actions          = Vector{InsertBeforeLoopHead}()
    insert_on_edge_actions                   = Vector{InsertOnEdge}()
    for action in actions
        if typeof(action) == InsertBeforeOrAfterStatement
            group_action(insert_before_or_after_statement_actions, action)
        elseif typeof(action) == InsertBeforeLoopHead
            group_action(insert_before_loop_head_actions, action)
        else
            group_action(insert_on_edge_actions, action)
        end
    end
    
    # Implement InsertOnEdge actions
    for action in insert_on_edge_actions
        if isempty(action.new_stmts)
            continue
        end

        new_bb, new_goto_stmt = CFGs.insertBetween(cfg, action.from_bb.label, action.to_bb.label)
        new_bb.statements     = action.new_stmts
        if new_goto_stmt != nothing
            push!(new_bb.statements, new_goto_stmt)
        end
    end

    # Implement InsertBeforeOrAfterStatement actions
    for action in insert_before_or_after_statement_actions
        if isempty(action.new_stmts)
            continue
        end

        for i = 1 : length(action.new_stmts)
            if action.before
                insert!(action.bb.statements, action.stmt_idx + i - 1, action.new_stmts[i])
            else
                insert!(action.bb.statements, action.stmt_idx + i, action.new_stmts[i])
            end
        end
    end

    # Implement InsertBeforeLoopHead actions
    for action in insert_before_loop_head_actions
        if isempty(action.new_stmts)
            continue
        end

        new_bb, new_goto_stmt = CFGs.insertBefore(cfg, action.loop.head, action.outside_loop, action.loop.back_edge)
        new_bb.statements     = action.new_stmts
        if new_goto_stmt != nothing
            push!(new_bb.statements, new_goto_stmt)
        end
    end
end
