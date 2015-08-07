# Distributivity analysis

@doc """ Is the expression an assignment? """
function is_assignment(expr :: Expr)
    head = expr.head
    args = expr.args
    return length(args) == 2 && (
            head == :(=) || head == :(+=) || head == :(-=) || head == :(*=) || 
            head == :(/=) || head == :(\=) || head == :(%=) || head == :(^=) ||
            head == :(&=) || head == :(|=) || head == :($=) || head == :(>>>=) ||
            head == :(>>=) || head == :(<<=)
           )
end

@doc """
Return true if the function call is distributive.
"""
function check_distributivity_of_function_call(
    ast         :: Expr,
    symbol_info :: Sym2TypeMap
)
    head = expr.head
    args = expr.args
        
    # A typical call are in the following forms
    #   Expr(:call, :*, :A, :x)
    #   Expr(:call, :(:call, top(getfield), SparseAccelerator,:SpMV), :A, :x)
    # The first argument is the function, the others are the arguments for it.
    arg_types = ntuple(i-> type_of_ast_node(args[i+1], symbol_info), length(args) - 1)
    all_numbers, some_arrays = are_numbers_or_arrays(expr.typ, arg_types)
    distributive = false
    if all_numbers || !some_arrays
        # The function call's result and arguments are all numbers, or 
        # some are non-numbers (like Range{UInt64}) but not regular arrays. 
        # No arrays to care about. 
        distributive = true
    else
        module_name, function_name = resolve_call_names(args)
        if function_name == ""
            throw(UnresolvedFunction(head, args[1]))
        end
        fd = lookup_function_description(module_name, function_name, arg_types)
        if fd != nothing
            distributive = fd.distributive
        else
            throw(UndescribedFunction(module_name, function_name, arg_types))
        end
    end
    
    if distributive
        # In case some args are trees themselves (even if the args' result types
        # are Number), we should check the args further
        for x in args[2 : end]
            if !check_distributivity(x, symbol_info)
                return false
            end
        end
        # ISSUE: what if an arg is an array element? It is a scalar, but it involves an array.
        return true
    else
        return false
    end
end

@doc """
Return true if reordering is distributive over the AST.
"""
function check_distributivity(
    ast         :: Any,
    symbol_info :: Sym2TypeMap
)
    local asttyp = typeof(ast)
    if asttyp <: Tuple
        for i = 1:length(ast)
            if !check_distributivity(ast[i], symbol_info)
                return false
            end
        end
        return true
    elseif asttyp == Expr
        local head = ast.head
        local args = ast.args
        if head == :lambda
            # TODO: LambdaHandling package should identify the body
            body = args[3]
            assert(typeof(body) == Expr)
            assert(typeof(body.head) == :body)
            return check_distributivity(body, symbol_info)
        elseif head == :body
            for i = 1:length(args)
                if !check_distributivity(args[i], symbol_info)
                    return false
                end
            end
            return true
        elseif is_assignment(ast)
            for i = 1:length(args)
                if !check_distributivity(args[i], symbol_info)
                    return false
                end
            end
            return true
        elseif head == :return
            return check_distributivity(args, symbol_info)
        elseif head == :call || head == :call1
            return check_distributivity_of_function_call(ast, symbol_info)
        elseif head == :gotoifnot
            if_clause  = args[1]
            return check_distributivity(if_clause, symbol_info)
        elseif head == :line
            return true
        else
            throw(UnknownExprDistributivity(head, args))
        end
    elseif asttyp == SymbolNode  || asttyp == Symbol   || asttyp == GenSym ||
           asttyp == LabelNode   || asttyp == GotoNode || asttyp == LineNumberNode ||
           asttyp == ASCIIString || asttyp == LambdaStaticData ||
           asttyp <: Number      || asttyp == NewvarNode
        return true
    elseif asttyp.name == Array.name 
        # There can be such arrays like [:(((top(getfield))(SparseAccelerator,:SpMV))(A,x)]
        # So we have to check each element of the array.
        for element in  ast
            if !check_distributivity(args[i], symbol_info)
                return false
            end
        end
        return true
    else
        throw(UnknownASTDistributivity(ast))
    end
end

@doc """
Return true if reordering is distributive over all operations in the region.
"""
function check_distributivity(
    region      :: LoopRegion,
    cfg         :: CFG, 
    symbol_info :: Sym2TypeMap
)
    try
        blocks  = cfg.basic_blocks
        for bb_index in region.loop.members
            bb = blocks[bb_index]
            for stmt in bb.statements
                distributive = check_distributivity(stmt.tls.expr, symbol_info)
                if !distributive
                    return false
                end
            end
        end
    catch ex
        dprintln(1, 0, "Exception during distributivity analysis! Loop region is skipped for reordering.")
        dprintln(1, 1, ex)
        return false
    end
    return true
end