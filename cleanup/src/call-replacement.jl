TypedExprNode(x...) = LivenessAnalysis.TypedExpr(x...)

@doc """ Check if two variables are the same. """
function same_variable(x :: Any, y :: Any)
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

@doc """ Check if a variable is an argument of an expression. """
function is_an_arg(x :: Any, ast :: Any)
    if typeof(ast) <: Expr
        for arg in ast.args
            if same_variable(x, arg)
                return true
            end
        end
    end
    return false
end

@doc """
Pre-processing function: check if the LHS of an assignment is in the RHS. 
"""
function LHS_in_RHS(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: String,
    fknob_deletor :: String
)
    assert(ast.head == :(=))
    return is_an_arg(ast.args[1], ast.args[2]);
end

# There are two kinds of patterns: one is ExprPattern, the other InPlaceUpdatePattern.
# ExprPattern is for matching an expression within a statement. InPlaceUpdatePattern
# is for matching two expressions in two adjacent statements.

@doc """
For an Expr AST node, describe it for pattern matching and replacement.
Skeleton is a tuple in the form of (head, arg1_type or arg1, arg2_type, arg3_type, ...).
The second element is either a type, or an argument, depending on head.
Examples skeletons: 
    (:call, GlobalRef(Main, :dot), Vector, Vector)
    (:=, Vector, Vector)

Each element (except head) in the skeleton might also be a pattern that can be 
pattern matched. Their patterns are described in sub_expr_patterns. This enables 
a match across several expressions. :NO_SUB_PATTERNS means there is no sub-expr
patterns to match. Otherwise, sub_expr_patterns is a tuple with the same size of
skeleton, and for each element, :nothing means no sub-expr-pattern for the Expr's
arg at that position, and any other pattern name means a sub pattern to match.

Substitute is a tuple showing how to replace the Expr according to the pattern.
Symbolic arguments are used here: :NO_CHANGE means that no replacement is needed
(the pattern is for matching only); :arg1 means the Expr's args[1], :aarg12 means
the Expr's args[1].args[2], :naarg12 means the negative :aarg12, etc.

Pre_ and post_processing are optional call backs before and after replace. Usually
they do nothing.

Now that the new AST (usually for a library function call) is generated, we need
to prepare he function-specific context info (fknob). Fknob_creator and _deletor
describe how to create and delete it. Usually they are empty, which means no
fknob. But if there is one, since the function may have matrix inputs/output,
the matrix-specific context info (mknobs) may need to be remembered in the fknob.
These matrices are specified by matrices_to_track. Symbolic arguments may be 
used here, in terms of the arguments in the new AST.
"""
immutable ExprPattern <: Pattern
    name              :: String # Name of the pattern
    skeleton          :: Tuple
    sub_expr_patterns :: Tuple
    pre_processing    :: Function
    substitute        :: Tuple
    post_processing   :: Function
    fknob_creator     :: String
    fknob_deletor     :: String
    matrices_to_track :: Tuple
end

@doc """
Match the following pattern that crosses two adjacent assignment statements: 
    y = f(...) # f is some function. y is a Symbol or GenSym
    result = y
Result must have already been allocated space. Otherwise, we cannot use the pattern.

We would replace the two statements into the following instead:
    f!(result, ...)
    y = result
ASSUMPTION: the only difference between f and f! is that f! has the result 
as its first parameter (Space already allocated), while f must allocate space
for its result before it returns. By changing f to f!, we save memory allocation.

This pattern is frequently seen from a source line like x += ...
"""
immutable InPlaceUpdatePattern <: Pattern
    name        :: String # Name of the pattern
    f_skeleton  :: Tuple  # Skeleton of f
    f!          :: Tuple  # The substitute
end

# Below are the expr_patterns we care about.

# Patterns that are used only for matching (sub-expressions).
const SpMV_3_parameters_pattern = ExprPattern(
    "SpMV_3_parameters_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     Number, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_4_parameters_pattern = ExprPattern(
    "SpMV_4_parameters_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     Number, SparseMatrixCSC, Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_6_parameters_pattern = ExprPattern(
    "SpMV_6_parameters_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const number_times_matrix_vector_pattern = ExprPattern(
    "number_times_matrix_vector_pattern",
    (:call, GlobalRef(Main, :*), Number, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE,),
    do_nothing,
    "",
    "",
    ()
)

const number_times_vector_pattern = ExprPattern(
    "number_times_vector_pattern",
    (:call, GlobalRef(Main, :*), Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const WAXPBY_4_parameters_pattern = ExprPattern(
    "WAXPBY_4_parameters_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     Number, Vector, Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

# Patterns that do transformation

const dot_pattern1 = ExprPattern(
    "dot_pattern1",
    (:call, GlobalRef(Main, :dot), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:dot)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_pattern1 = ExprPattern(
    "SpMV_pattern1",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_pattern2 = ExprPattern(
    "SpMV_pattern2",
    (:call, GlobalRef(Main, :A_mul_B!), Number, SparseMatrixCSC, Vector, Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg6, :arg2, :arg3, :arg4, :arg5, :arg6, 0.0),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_pattern3 = ExprPattern(
    "SpMV_pattern3",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, number_times_matrix_vector_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_pattern4 = ExprPattern(
    "SpMV_pattern4",
    (:call, GlobalRef(Main, :*), Number, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3, :arg4),
    do_nothing,
    "",
    "",
    ()
)

const SpMV_pattern5 = ExprPattern(
    "SpMV_pattern5",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, SpMV_3_parameters_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3),
    do_nothing,
    "",
    "",
    ()
)

const SpMV!_pattern1 = ExprPattern(
    "SpMV!_pattern1",
    (:call, GlobalRef(Main, :A_mul_B!), Vector, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg2, :arg3, :arg4),
    do_nothing,
    "",
    "",
    ()
)

const SpMV!_pattern2 = ExprPattern(
    "SpMV!_pattern2",
    (:(=), Vector, Vector),
    (nothing, nothing, SpMV_6_parameters_pattern),
    LHS_in_RHS,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg1, :aarg22, :aarg23, :aarg24, :aarg25, :aarg26, :aarg27),
    do_nothing,
    "",
    "",
    ()
)

const SpMV!_pattern3 = ExprPattern(
    "SpMV!_pattern3",
    (:(=), Vector, Vector),
    (nothing, nothing, SpMV_3_parameters_pattern),
    LHS_in_RHS,
    (TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg1, :aarg22, :aarg23, :aarg24, 0, :arg1, :aarg27),
    do_nothing,
    "",
    "",
    ()
)

const WAXPBY_pattern1 = ExprPattern(
    "WAXPBY_pattern1",
    (:call, GlobalRef(Main, :+), Vector, Vector),
    (nothing, nothing, nothing, number_times_vector_pattern),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     1, :arg2, :aarg32, :aarg33),
    do_nothing,
    "",
    "",
    ()
)

const WAXPBY_pattern2 = ExprPattern(
    "WAXPBY_pattern2",
    (:call, GlobalRef(Main, :-), Vector, Vector),
    (nothing, nothing, nothing, number_times_vector_pattern),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     1, :arg2, :naarg32, :aarg33),
    do_nothing,
    "",
    "",
    ()
)

const WAXPBY!_pattern1 = ExprPattern(
    "WAXPBY!_pattern1",
    (:(=), Vector, Vector),
    (nothing, nothing, WAXPBY_4_parameters_pattern),
    LHS_in_RHS,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY!)),
     :arg1, :aarg22, :aarg23, :aarg24, :aarg25),
    do_nothing,
    "",
    "",
    ()
)

const element_wise_multiply_pattern1 = ExprPattern(
    "element_wise_multiply_pattern1",
    (:call, GlobalRef(Main, :.*), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:element_wise_multiply)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    ()
)

expr_patterns = [
    dot_pattern1,
    #WAXPBY_pattern,
    SpMV_pattern1,
    SpMV_pattern2,
    SpMV_pattern3,
    SpMV_pattern4,
    SpMV_pattern5,
    SpMV!_pattern1,
    SpMV!_pattern2,
    WAXPBY_pattern1,
    WAXPBY_pattern2,
    WAXPBY!_pattern1,
    element_wise_multiply_pattern1
    #SpMV_pattern2,
    #WAXPBY!_pattern,
]

const SpMV!_in_place_update_pattern1 = InPlaceUpdatePattern(
    "SpMV!_in_place_update_pattern1",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :result, :arg2, :arg3, :arg4, :arg5, :arg6, :arg7)
)

const WAXPBY!_in_place_update_pattern1 = InPlaceUpdatePattern(
    "WAXPBY!_in_place_update_pattern1",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
     Number, Vector, Number, Vector),
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY!)),
     :result, :arg2, :arg3, :arg4, :arg5)
)

in_place_update_patterns = [
    SpMV!_in_place_update_pattern1,
    WAXPBY!_in_place_update_pattern1
]

@doc """ Return an Expr AST node's skeleton """
function expr_skeleton(
    ast         :: Expr, 
    symbol_info :: Sym2TypeMap
)
    head = ast.head
    # So far, we handle only call and assignment expressions. But there might be
    # other kinds of expressions like LineNumber. For them, make the second
    # parameter the real arg: it is not matched anyway.
    args = ast.args
    skeleton =  ntuple(i-> (
                             (i == 1) ? head
                                      : (i == 2) ? (head == :(=) ? type_of_ast_node(args[1], symbol_info) : args[1])
                                                 : type_of_ast_node(args[i - 1], symbol_info)
                           ), length(args) + 1)
    skeleton
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
            arg = TypedExprNode(type_of_ast_node(arg, symbol_info), :call, :(-), arg)
        end
    end
    arg
end

replacement_arg(arg :: Any, args :: Vector, symbol_info :: Sym2TypeMap) =
    replacement_arg(arg, args, nothing, symbol_info)

@doc """ Replace the AST based on the expression pattern. """
function replace(
    pattern     :: ExprPattern,
    ast         :: Expr,
    symbol_info :: Sym2TypeMap
)
    substitute = pattern.substitute
    if length(substitute) == 1 && substitute[1] == :NO_CHANGE
        return
    end

    orig_ast = copy(ast)
    ast.head = substitute[1]
    empty!(ast.args)
    for i = 2 : length(substitute)
        arg = replacement_arg(substitute[i], orig_ast.args, symbol_info)
        push!(ast.args, arg)
    end
end

@doc """ 
Match two skeletons. Return true if matched. 
"""
function match_skeletons(
    skeleton1 :: Tuple,
    skeleton2 :: Tuple
)
    if length(skeleton1) == length(skeleton2)
        if (skeleton1[1] == skeleton2[1]) && 
           ((skeleton1[1] == :call && skeleton1[2] == skeleton2[2]) ||
            (skeleton1[1] == :(=)  && skeleton1[2] <: skeleton2[2]))
            for i in 3:length(skeleton1)
                if !(skeleton1[i] <: skeleton2[i])
                    return false
                end
            end
            return true
        end
    end
    return false
end

@doc """ 
Match and replace a pattern. Return true if matched. 
Note: match and replace has to be done in the same time, because one expression
may need its sub-expressions be matched and replaced first.
"""
function match_replace(
    pattern    :: ExprPattern,
    skeleton   :: Tuple,
    ast        :: Expr,
    call_sites :: CallSites
)
    #if pattern == CS_last_resort_pattern
    #    # This is the only pattern we do special handling. When it is reached, 
    #    # it means all other patterns have been tried and do not work, and this
    #    # one will tell us to do something special
    #    return pattern.pre_processing(ast, call_sites, pattern.fknob_creator, pattern.fknob_deletor)
    #end
    
    symbol_info = call_sites.symbol_info
    if match_skeletons(skeleton, pattern.skeleton)
        # Check sub-expr_patterns
        if length(pattern.sub_expr_patterns) == 1 && 
           pattern.sub_expr_patterns[1] == :NO_SUB_PATTERNS
           # Do nothing
        else
            assert(pattern.sub_expr_patterns[1] == nothing)
            for i = 2 : length(pattern.sub_expr_patterns)
                sub_pattern = pattern.sub_expr_patterns[i]
                if sub_pattern != nothing
                    if typeof(ast.args[i- 1]) <: Expr
                        sub_expr_skeleton = expr_skeleton(ast.args[i- 1], symbol_info)
                        if !match_replace(sub_pattern, sub_expr_skeleton, ast.args[i - 1], call_sites)
                            return false
                        end
                    else
                        return false
                    end
                end
            end
        end
        
        if pattern.pre_processing != do_nothing
            if !pattern.pre_processing(ast, call_sites, pattern.fknob_creator,
                                       pattern.fknob_deletor)
                return false
            end
        end
        
        dprintln(1, 0, "\n\nReplace")
        dsprintln(1, 1, symbol_info, ast)

        replace(pattern, ast, symbol_info)
        
        if pattern.post_processing != do_nothing
            if !pattern.post_processing(ast, call_sites, pattern.fknob_creator,
                                        pattern.fknob_deletor, pattern.matrices_to_track)
                # AST has already been changed by replace(). However, post 
                # processing fails. That AST might be wrong. So abort 
                dprintln(1, 0, "to")
                dsprintln(1, 1, symbol_info, ast)
                dprintln(1, 0, "according to pattern but post-processing failed.")
                dprintln(1, 1, pattern)
                throw(PostPatternReplacementFailure(pattern))
            end
        end
        
        dprintln(1, 0, "to")
        dsprintln(1, 1, symbol_info, ast)
        dprintln(1, 0, "according to pattern")
        dprintln(1, 1, pattern)
        
        return true
    end
    return false
end

@doc """
Match an expression pattern and do replacement.
"""
function match_replace_an_expr_pattern(
    ast        :: Any,
    call_sites :: CallSites
)
    symbol_info = call_sites.symbol_info
    patterns    = call_sites.patterns
    if typeof(ast) <: Expr
        # Match against each arg first. Replacement may happen on the args.
        for arg in ast.args
            match_replace_an_expr_pattern(arg, call_sites)
        end

        # Now match against the Expr. Replacement may happen on the expression.
        skeleton = expr_skeleton(ast, symbol_info)
        for pattern in patterns
            success = match_replace(pattern, skeleton, ast, call_sites)
            if success
                return nothing
            end
        end
    end
    return nothing
end

@doc """
Match an expression pattern and do replacement.
"""
function match_replace_an_expr_pattern(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    match_replace_an_expr_pattern(ast, call_sites)
end

@doc """ Replace two Expr nodes based on the in-place-update pattern. """
function replace(
    pattern     :: InPlaceUpdatePattern,
    prev_expr   :: Expr,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap
)
    dprintln(1, 0, "Replace")
    dsprintln(1, 1, symbol_info, prev_expr)
    dsprintln(1, 1, symbol_info, expr)

    # Replace prev_expr with the f!
    result         = copy(expr.args[1])
    f              = copy(prev_expr.args[2].args)
    f!             = pattern.f!
    prev_expr.head = f![1]
    empty!(prev_expr.args)
    for i = 2 : length(f!)
        arg = replacement_arg(f![i], f, result, symbol_info)
        push!(prev_expr.args, arg)
    end

    # Swith the two arguments of expr.
    arg = copy(expr.args[1])
    expr.args[1] = expr.args[2]
    expr.args[2] = arg
    
    dprintln(1, 0, "to")
    dsprintln(1, 1, symbol_info, prev_expr)
    dsprintln(1, 1, symbol_info, expr)
    dprintln(1, 0, "according to pattern")
    dprintln(1, 1, pattern)
end

@doc """
Check if the two expression, from adjacent statements, are of the following pattern:
    y = f(...) # f is some function. y is a Symbol or GenSym
    result = y
Result must have already been allocated space.

Replace it into the following instead:
    f!(result, ...)
    y = result
"""
function match_replace_an_in_place_update_pattern(
    prev_expr   :: Expr,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap
)
    if prev_expr.head == :(=) && expr.head == :(=) &&
       length(prev_expr.args) == 2 && length(expr.args) == 2 &&
       prev_expr.args[1] == expr.args[2] && 
       (typeof(expr.args[2]) == Symbol || typeof(expr.args[2]) == GenSym) &&
       typeof(prev_expr.args[2]) == Expr
       
        # Check that the result is one of the input parameter of f, in which case
        # it must have been allocated space already.
        if !is_an_arg(expr.args[1], prev_expr.args[2])
            return false
        end

        f_skeleton = expr_skeleton(prev_expr.args[2], symbol_info)
        for pattern in in_place_update_patterns
            if match_skeletons(f_skeleton, pattern.f_skeleton)
                replace(pattern, prev_expr, expr, symbol_info)
                return true
            end
        end
    end
    return false
end

@doc """ 
Pattern match and replace the code that is functionally equivalent to SpMV, SpMV!,
dot, WAXPBY, WAXPBY!, etc. with calls to the corresponding SPMP library functions.
"""
function replace_calls(
    symbol_info :: Sym2TypeMap, 
    cfg         :: CFG
)
    call_sites  = CallSites(Set{CallSite}(), WholeFunction(), symbol_info, 
                            expr_patterns, Vector{Action}(), nothing)
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
                # Between this and the previous expression, try to optimize for an in-place update
                match_replace_an_in_place_update_pattern(prev_expr, expr, symbol_info)
            end
            
            prev_expr = expr
        end
    end
end
