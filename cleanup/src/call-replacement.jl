abstract Pattern

@doc """
For an Expr AST node, describe it for pattern matching and replacement.
"""
immutable ExprPattern <: Pattern
    head              :: Symbol # A function call in general
    signature         :: Tuple  # Tuple{Symbol or Type}. Signature of the function call
    sub_expr_patterns :: Tuple  # Tuple{ExprPattern}. Each arg in the signature might be an sub-expression that can be pattern matched. This enables a match across several expressions
    new_args          :: Tuple  # Tuple{Any}. The new args to replace the matched args
end

@doc """
Match the following pattern that crosses two adjacent statements: 
    y = f_x   # f_x is some function with x as an input. y is a Symbol or GenSym
    x = y
We would replace the two statements into the following instead:
    new_f_x
    y = x
new_f_x is equivalent to x = f_x, but is an in_place update function, thus 
saving memory allocation for the temporary result of f_x.

This pattern is frequently seen from a source line like x += ...
"""
immutable SymPattern <: Pattern
    f_x_signature :: Tuple  # Tuple{Symbol or Type}. Signature of f_x
    x_idx         :: Int    # Position of x in f_x_signature
    new_f_x       :: Tuple  # Tuple{Any}. The args of the new expression
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

const SpMV!_pattern1 = ExprPattern(
    :call,
    (GlobalRef(Main, :A_mul_B!), Vector, SparseMatrixCSC, Vector),
    (nothing, nothing, nothing),
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg2, :arg3, :arg4)
)

const number_times_vector_pattern = ExprPattern(
    :call,
    (GlobalRef(Main, :*), Number, Vector),
    (nothing, nothing, nothing),
    (:arg1, :arg2, :arg3)
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

expr_patterns = [
    dot_pattern,
    #WAXPBY_pattern,
    SpMV_pattern1,
    SpMV_pattern2,
    SpMV!_pattern1,
    WAXPBY_pattern1,
    WAXPBY_pattern2
    #SpMV_pattern2,
    #WAXPBY!_pattern,
]

const WAXPBY!_pattern = SymPattern(
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     Number, Vector, Number, Vector),
    3,
    (LivenessAnalysis.TypedExpr(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY!)),
     :arg3, :arg2, :arg3, :arg4, :arg5)
)

sym_patterns = [
    WAXPBY!_pattern
]

@doc """ Return an Expr AST node's head and signature """
function expr_head_signature(
    ast         :: Expr, 
    symbol_info :: Sym2TypeMap
)
    head = ast.head
    args = ast.args
    signature =  ntuple(i-> (i == 1 ? args[1] : 
        type_of_ast_node(args[i], symbol_info)), length(args))
    head, signature
end

@doc """ 
The real arg to replace the symbolic arg (like :arg1) in the patterns.
"""
function replacement_arg(
    arg         :: Any, 
    args        :: Vector, 
    symbol_info :: Sym2TypeMap
)
    if typeof(arg) == Symbol
        arg_string = string(arg)
        if length(arg_string) > 3 && arg_string[1 : 3] == "arg" 
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

@doc """ Replace the AST based on the pattern. """
function replace(
    pattern     :: ExprPattern,
    ast         :: Expr,
    symbol_info :: Sym2TypeMap
)
    new_args = pattern.new_args
    args     = copy(ast.args)
    for i = 1 : length(new_args)
        arg = replacement_arg(new_args[i], args, symbol_info)
        if i <= length(ast.args)
            ast.args[i] = arg
        else
            push!(ast.args, arg)
        end
    end
    
    println("after replace new ast=", ast)
end

@doc """ 
Match and replace a pattern. Return true if matched. 
Note: match and replace must be do in the same time, because one expression may
need its sub-expressions be matched and replaced first.
"""
function match_replace(
    pattern     :: ExprPattern,
    head        :: Symbol,
    signature   :: Tuple,
    ast         :: Expr,
    symbol_info :: Sym2TypeMap
)
if pattern == WAXPBY_pattern1 || pattern == number_times_vector_pattern
println("pattern is=", pattern)
println("\tto match =", head, " ", signature, "\n\tast=", ast)

end

    if pattern.head == head && length(pattern.signature) == length(signature)
        if signature[1] == pattern.signature[1]
            for i in 2:length(signature) 
                if !(signature[i] <: pattern.signature[i])
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
        println(".... success match")
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
        head, signature = expr_head_signature(ast, call_sites.symbol_info)
        for pattern in expr_patterns
            success = match_replace(pattern, head, signature, ast, call_sites.symbol_info)
            if success
                return nothing
            end
        end
    end
    return nothing
end

@doc """ Replace the two Expr nodes based on the pattern. """
function replace(
    pattern     :: SymPattern,
    prev_expr   :: Expr,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap
)
    # Replace prev_expr with the new_f_x
    f_x     = copy(prev_expr.args[2].args)
    new_f_x = pattern.new_f_x
    for i = 1 : length(new_f_x)
        arg = replacement_arg(new_f_x[i], f_x, symbol_info)
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

@doc """ Check if two variables are the same. """
function same_symbol(x :: Any, y :: Any)
    x1 = x
    y1 = y
    if typeof(x) <: SymbolNode
        x1 = x.name
    end
    if typeof(y) <: SymbolNode
        y1 = y.name
    end
    x1 == y1
end
 
@doc """
Check if the two expression, from adjacent statements, are of the following pattern:
    y = f_x   # f_x is some function with x as an input. y is a Symbol or GenSym
    x = y
Replace it into the following instead:
    new_f_x
    y = x
"""
function match_replace_an_sym_pattern(
    prev_expr   :: Expr,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap
)
println("mr sympatterh: prev=", prev_expr)
println("\t expr=", expr)
    if prev_expr.head == :(=) && expr.head == :(=) &&
       prev_expr.args[1] == expr.args[2] && 
       length(prev_expr.args) == 2 && length(expr.args) == 2 &&
       typeof(prev_expr.args[2]) == Expr &&
       (typeof(expr.args[2]) == Symbol || typeof(expr.args[2]) == GenSym)
       
       println("prev2.head=", prev_expr.args[2].head)
       println("prev2.args=", prev_expr.args[2].args)
       
        for pattern in sym_patterns
        println("\t\ttry pattern=", pattern)
            f_x       = prev_expr.args[2].args
println("f_x is:", f_x)
println("\tf_x type is:", typeof(f_x))
            x         = expr.args[1]

println("****:",  length(f_x), "  ", length(pattern.f_x_signature))
println("^^^^" , x,  "typeofx=", typeof(x), "  ",  f_x[pattern.x_idx],  "typeofFx=", typeof(f_x[pattern.x_idx]))
println("%%%", x == f_x[pattern.x_idx])

            if length(f_x) == length(pattern.f_x_signature) && same_symbol(x, f_x[pattern.x_idx])
                if f_x[1] == pattern.f_x_signature[1]
                    for i in 2:length(pattern.f_x_signature) 
                    println("~~~ ", i, " ", "f_x[i]=", f_x[i], " type=",  typeof(f_x[i]) )
                    println("\t~~~ sigf_x[i]=", pattern.f_x_signature[i], " type=",  typeof(pattern.f_x_signature[i]) )
                    
                        if !(type_of_ast_node(f_x[i], symbol_info) <: pattern.f_x_signature[i])
                println("\t\t compare ", i, " failed")
                            return false
                        end
                    end
                else
                println("\t\t compare 1 failed")
                    return false
                end
                
                println(".... Sympatternsuccess match")
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
                match_replace_an_sym_pattern(prev_expr, expr, symbol_info)
            end
            
            prev_expr = expr
        end
    end
end