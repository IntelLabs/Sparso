@doc """
For an Expr AST node, describe it for pattern matching and replacement.
"""
immutable Pattern
    head            :: Symbol # A function call in general
    signature       :: Tuple  # Signature of the function call
    substitute_args :: Tuple  # The new args to replace the matched args
end

# Below are the patterns we care about. Here "arg2", etc. represents the original
# args[2], etc.
const dot_pattern = Pattern(
    :call,
    (GlobalRef(Main, :dot), Vector, Vector),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:Dot)),
     :arg2, :arg3)
)

patterns = [
    dot_pattern
    #,
    #WAXPBY_pattern,
    #SpMV_pattern1,
    #SpMV_pattern2,
    #WAXPBY!_pattern,
]

@doc """ Match a pattern. Return true if matched. """
function match(
    pattern   :: Pattern,
    head      :: Symbol,
    signature :: Tuple
)
    if pattern.head == head && length(pattern.signature) == length(signature)
        matched = true
        if signature[1] == pattern.signature[1]
            for i in 2:length(signature) 
                if !(signature[i] <: pattern.signature[i])
                    matched = false
                    break
                end
            end
        else
            matched = false
        end
        return matched
    end
    return false
end

@doc """ Replace the AST based on the pattern. """
function replace(
    pattern :: Pattern,
    ast     :: Expr
)
    substitute_args = pattern.substitute_args
    args            = copy(ast.args)
    for i = 1 : length(substitute_args)
        arg = substitute_args[i]
        if string(arg)[1 : 3] == "arg" 
            # Replace the special symbols like :arg2 with real arguments.
            if arg == :arg2
                arg = args[2]
            elseif arg == :arg3
                arg = args[3]
            elseif arg == :arg4
                arg = args[4]
            elseif arg == :arg5
                arg = args[5]
            else
                throw(UnhandledSubstituteArg(arg, pattern))
            end
        end
    
        if i <= length(ast.args)
            ast.args[i] = arg
        else
            push!(ast.args, arg)
        end
    end
end

@doc """
Match a pattern and do replacement.
"""
function match_a_pattern(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    if typeof(ast) <: Expr
        head = ast.head
        args = ast.args
        signature =  ntuple(i-> (i == 1 ? args[1] : 
            type_of_ast_node(args[i], call_sites.symbol_info)), length(args))
        for pattern in patterns
            success = match(pattern, head, signature)
            if success
                replace(pattern, ast)
                return nothing
            end
        end
    end
    return nothing
end

@doc """ 
Pattern match and replace the code that is functionally equivalent to SpMV, dot,
WAXPBY, etc. with calls to the corresponding SPMP library functions.
"""
function replace_calls(
    symbol_info :: Sym2TypeMap, 
    cfg         :: CFG
)
    call_sites  = CallSites(Set{CallSite}(), symbol_info)
    for (bb_idx, bb) in cfg.basic_blocks
        for stmt in bb.statements
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            CompilerTools.AstWalker.AstWalk(expr, match_a_pattern, call_sites)
        end
    end
end