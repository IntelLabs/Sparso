TypedExprNode(x...) = LivenessAnalysis.TypedExpr(x...)

abstract Pattern

@doc """
For an Expr AST node, describe it for pattern matching and replacement.
"""
immutable ExprPattern <: Pattern
    name              :: String # Name of the pattern
    head              :: Symbol # A function call in general
    skeleton         :: Tuple  # Tuple{Symbol or Type}. Signature of the function call
    sub_expr_patterns :: Tuple  # Tuple{ExprPattern}. Each arg in the skeleton might be an sub-expression that can be pattern matched. This enables a match across several expressions
    new_args          :: Tuple  # Tuple{Any}. The new args to replace the matched args
end

@doc """
Match the following pattern that crosses two adjacent statements: 
    y = f(...) # f is some function. y is a Symbol or GenSym
    result = y
We would replace the two statements into the following instead:
    result = f!(result, ...)
    y = result
ASSUMPTION: the only difference between f and f! is that f! has the result 
as its first parameter (Space already allocated), while f must allocate space
for its result before it returns. By changing f to f!, we save memory allocation.

This pattern is frequently seen from a source line like x += ...
"""
immutable InPlaceUpdatePattern <: Pattern
    f_signature :: Tuple  # Tuple{Symbol or Type}. Signature of f
    f!          :: Tuple  # Tuple{Any}.
end

# Below are the expr_patterns we care about. Here "arg2", etc. represents the original
# args[2], etc.
const dot_pattern = ExprPattern(
    :call,
    (GlobalRef(Main, :dot), Vector, Vector),
    (nothing, nothing, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:Dot)),
     :arg2, :arg3)
)

const SpMV_pattern1 = ExprPattern(
    :call,
    (GlobalRef(Main, :*), SparseMatrixCSC, Vector),
    (nothing, nothing, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3)
)

const SpMV_pattern2 = ExprPattern(
    :call,
    (GlobalRef(Main, :A_mul_B!), Number, SparseMatrixCSC, Vector, Number, Vector),
    (nothing, nothing, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3, :arg4, :arg5, :arg6)
)

const number_times_matrix_vector_pattern = ExprPattern(
    :call,
    (GlobalRef(Main, :*), Number, SparseMatrixCSC, Vector),
    (nothing, nothing, nothing),
    (:arg1, :arg2, :arg3)
)

const SpMV_pattern3 = ExprPattern(
    :call,
    (GlobalRef(Main, :+), Vector, Number),
    (nothing, number_times_matrix_vector_pattern, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3)
)

const SpMV!_pattern1 = ExprPattern(
    :call,
    (GlobalRef(Main, :A_mul_B!), Vector, SparseMatrixCSC, Vector),
    (nothing, nothing, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg2, :arg3, :arg4)
)

const SpMV_6_parameters_pattern = ExprPattern(
    :call,
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    (nothing, nothing, nothing, nothing, nothing, nothing),
    (:NO_CHANGE, )
)

const SpMV!_pattern2 = ExprPattern(
    :(=),
    (Vector, Vector),
    (nothing, SpMV_6_parameters_pattern),
    (TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg1, :aarg22, :aarg23, :aarg24, :aarg25, :aarg26, :aarg27)
)

const number_times_vector_pattern = ExprPattern(
    :call,
    (GlobalRef(Main, :*), Number, Vector),
    (nothing, nothing, nothing),
    (:NO_CHANGE, )
)

const WAXPBY_pattern1 = ExprPattern(
    :call,
    (GlobalRef(Main, :+), Vector, Vector),
    (nothing, nothing, number_times_vector_pattern),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     1, :arg2, :aarg32, :aarg33)
)

const WAXPBY_pattern2 = ExprPattern(
    :call,
    (GlobalRef(Main, :-), Vector, Vector),
    (nothing, nothing, number_times_vector_pattern),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     1, :arg2, :naarg32, :aarg33)
)

const PointwiseMultiply_pattern1 = ExprPattern(
    :call,
    (GlobalRef(Main, :.*), Vector, Vector),
    (nothing, nothing, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:PointwiseMultiply)),
     :arg2, :arg3)
)


expr_patterns = [
    dot_pattern,
    #WAXPBY_pattern,
    SpMV_pattern1,
    SpMV_pattern2,
    SpMV_pattern3,
    SpMV!_pattern1,
    SpMV!_pattern2,
    WAXPBY_pattern1,
    WAXPBY_pattern2,
    PointwiseMultiply_pattern1
    #SpMV_pattern2,
    #WAXPBY!_pattern,
]

const SpMV!_in_place_update_pattern1 = InPlaceUpdatePattern(
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :result, :arg2, :arg3, :arg4, :arg5, :arg6, :arg7)
)

const WAXPBY!_in_place_update_pattern1 = InPlaceUpdatePattern(
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     Number, Vector, Number, Vector),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY!)),
     :result, :arg2, :arg3, :arg4, :arg5)
)

in_place_update_patterns = [
    SpMV!_in_place_update_pattern1,
    WAXPBY!_in_place_update_pattern1
]

@doc """ Return an Expr AST node's head and skeleton """
function expr_head_signature(
    ast         :: Expr, 
    symbol_info :: Sym2TypeMap
)
    head = ast.head
    args = ast.args
    skeleton =  ntuple(i-> (i == 1 ? args[1] : 
        type_of_ast_node(args[i], symbol_info)), length(args))
    head, skeleton
end

@doc """ 
The real arg to replace the symbolic arg (like :arg1) in the patterns.
"""
function replacement_arg(
    arg         :: Any, 
    args        :: Vector,
    result      :: Any,       # For in_place update pattern only
    symbol_info :: Sym2TypeMap
)
    if typeof(arg) == Symbol
        arg_string = string(arg)
        if length(arg_string) == 6 && arg_string[1 : 6] == "result" 
            # This refers to the result arg of expr in an in_place update pattern.
            assert(result != nothing)
            arg = result
        elseif length(arg_string) > 3 && arg_string[1 : 3] == "arg" 
            # This refers to an arg in the form of argx. For
            # example, arg3 means args[3]
            x   = Int(arg_string[4]) - Int('0')
            arg = args[x]
        elseif length(arg_string) > 4 && arg_string[1 : 4] == "aarg"
            # This refers to a sub expression's arg in the form of aargxy. For
            # example, aarg31 means args[3].args[1]
            x = Int(arg_string[5]) - Int('0')
            y = Int(arg_string[6]) - Int('0')
            arg = args[x].args[y]
        elseif length(arg_string) > 5 && arg_string[1 : 5] == "naarg"
            # This refers to a sub expression's arg in the form of aargxy. For
            # example, aarg31 means args[3].args[1]
            x   = Int(arg_string[6]) - Int('0')
            y   = Int(arg_string[7]) - Int('0')
            arg = args[x].args[y]
            arg = LivenessAnalysis.TypedExpr(type_of_ast_node(arg, symbol_info), :call, :(-), arg)
        end
    end
    arg
end

replacement_arg(arg :: Any, args :: Vector, symbol_info :: Sym2TypeMap) =
    replacement_arg(arg, args, nothing, symbol_info)

@doc """ Replace the AST based on the pattern. """
function replace(
    pattern     :: ExprPattern,
    ast         :: Expr,
    symbol_info :: Sym2TypeMap
)
    new_args = pattern.new_args
    if length(new_args) == 1 && new_args[1] == :NO_CHANGE
        return
    end

    dprintln(1, 0, "Replace")
    dsprintln(1, 1, ast)

    args = copy(ast.args)
    for i = 1 : length(new_args)
        arg = replacement_arg(new_args[i], args, symbol_info)
        if i <= length(ast.args)
            ast.args[i] = arg
        else
            push!(ast.args, arg)
        end
    end
    
    dprintln(1, 0, "to")
    dsprintln(1, 1, ast)
    dprintln(1, 0, "according to pattern")
    dprintln(1, 1, pattern)
end

@doc """ 
Match and replace a pattern. Return true if matched. 
Note: match and replace must be do in the same time, because one expression may
need its sub-expressions be matched and replaced first.
"""
function match_replace(
    pattern     :: ExprPattern,
    head        :: Symbol,
    skeleton   :: Tuple,
    ast         :: Expr,
    symbol_info :: Sym2TypeMap
)
if pattern.head == :(=)  
println("matching ", ast,  " against ")
println("\t ", pattern)
println(" \t\ttype1=", type_of_ast_node(skeleton[1], symbol_info), "  ",  pattern.skeleton[1])

    if type_of_ast_node(skeleton[1], symbol_info) <: pattern.skeleton[1]
        println("\t\tsucc!")
    else 
        println("\t\tfail!")
    end
end
    if pattern.head == head && length(pattern.skeleton) == length(skeleton)
        if (pattern.head == :call && skeleton[1] == pattern.skeleton[1]) ||
           (pattern.head == :(=)  && type_of_ast_node(skeleton[1], symbol_info) <: pattern.skeleton[1])
            for i in 2:length(skeleton) 
                if !(skeleton[i] <: pattern.skeleton[i])
                    return false
                end
            end
        else
            return false
        end
        
        # Check sub-expr_patterns
        for i = 1 : length(pattern.sub_expr_patterns)
            sub_pattern = pattern.sub_expr_patterns[i]
            if sub_pattern != nothing
                sub_expr_head, sub_expr_signature = expr_head_signature(ast.args[i], symbol_info)
                if !match_replace(sub_pattern, sub_expr_head, sub_expr_signature, ast.args[i], symbol_info)
                    return false
                end
            end
        end
        replace(pattern, ast, symbol_info)
        return true
    end
    return false
end

@doc """
Match a pattern and do replacement.
"""
function match_replace_an_expr_pattern(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    if typeof(ast) <: Expr
        # Match against each arg first.
        for arg in ast.args
            match_replace_an_expr_pattern(arg, call_sites, top_level_number, is_top_level, read)
        end
        
        # Now match against the Expr.
        head, skeleton = expr_head_signature(ast, call_sites.symbol_info)
        for pattern in expr_patterns
            success = match_replace(pattern, head, skeleton, ast, call_sites.symbol_info)
            if success
                return nothing
            end
        end
    end
    return nothing
end

@doc """ Replace the two Expr nodes based on the pattern. """
function replace(
    pattern     :: InPlaceUpdatePattern,
    prev_expr   :: Expr,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap
)
    # Replace prev_expr with the f!
    result  = copy(expr.args[1])
    f       = copy(prev_expr.args[2].args)
    f!      = pattern.f!
    for i = 1 : length(f!)
        arg = replacement_arg(f![i], f, result, symbol_info)
        if i <= length(prev_expr.args[2].args)
            prev_expr.args[2].args[i] = arg
        else
            push!(prev_expr.args[2].args, arg)
        end
    end
    prev_expr.args[1] = expr.args[1]

    # Swith the two arguments of expr.
    arg = copy(expr.args[1])
    expr.args[1] = expr.args[2]
    expr.args[2] = arg
    
    println("Sympatern replace new prev_expr=", prev_expr)
    println("Sympatern replace new expr=", expr)
end

@doc """
Check if the two expression, from adjacent statements, are of the following pattern:
    y = f(...) # f is some function. y is a Symbol or GenSym
    result = y
Replace it into the following instead:
    result = f!(result, ...)
    y = result
"""
function match_replace_an_in_place_update_pattern(
    prev_expr   :: Expr,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap
)
    if prev_expr.head == :(=) && expr.head == :(=) &&
       prev_expr.args[1] == expr.args[2] && 
       (typeof(expr.args[2]) == Symbol || typeof(expr.args[2]) == GenSym) &&
       length(prev_expr.args) == 2 && length(expr.args) == 2 &&
       typeof(prev_expr.args[2]) == Expr
       
        f = prev_expr.args[2].args
        for pattern in in_place_update_patterns
            if length(f) == length(pattern.f_signature)
                if f[1] == pattern.f_signature[1]
                    for i in 2:length(pattern.f_signature) 
                        if !(type_of_ast_node(f[i], symbol_info) <: pattern.f_signature[i])
                            return false
                        end
                    end
                else
                    return false
                end
                replace(pattern, prev_expr, expr, symbol_info)
                return true
            end
            return false
        end
    end
    return false
end

@doc """ 
Pattern match and replace the code that is functionally equivalent to SpMV, dot,
WAXPBY, etc. with calls to the corresponding SPMP library functions. Also optimize
for GenSym expr_patterns.
"""
function replace_calls(
    symbol_info :: Sym2TypeMap, 
    cfg         :: CFG
)
    call_sites  = CallSites(Set{CallSite}(), symbol_info)
    for (bb_idx, bb) in cfg.basic_blocks
        prev_expr = nothing
        for stmt in bb.statements
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end
            
            # Try to pattern match and replace this expression with ExprPatterns.
            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
            
            if prev_expr != nothing
                # Between this and the previous expression, try to optimize for SymPatterns
                match_replace_an_in_place_update_pattern(prev_expr, expr, symbol_info)
            end
            
            prev_expr = expr
        end
    end
end