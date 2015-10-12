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
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString
)
    assert(ast.head == :(=))
    return is_an_arg(ast.args[1], ast.args[2]);
end

# There are two kinds of patterns: one is ExprPattern, the other TwoStatementsPattern.
# ExprPattern is for matching an expression within a statement (It can be the
# whole or part of the statement). TwoStatementsPattern is for matching
# two whole adjacent statements.

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
The matched function call can decide reordering or not, based on its power. If
it decides to reorder, the arrays it first reorders (First Arrays Reordered) are
specified by reoredering_FAR; the first of the arrays is the one that determines
permutation vectors, called Seed.
"""
immutable ExprPattern <: Pattern
    name              :: AbstractString # Name of the pattern
    skeleton          :: Tuple
    sub_expr_patterns :: Tuple
    pre_processing    :: Function
    substitute        :: Tuple
    post_processing   :: Function
    fknob_creator     :: AbstractString
    fknob_deletor     :: AbstractString
    matrices_to_track :: Tuple
    reordering_power  :: Int
    reordering_FAR    :: Tuple
end

@doc """
Match a pattern that crosses two adjacent statements. The fields are
the same as ExprPattern, except that it has two skeletons to match, and 
two new skeletons for replacement. See the description of ExprPattern for the
common fields' meaning. 

After replacement, the new first statement is the one that has the total
effect of the original two statements, while the new second statement is
usually just a copy (maybe useless). So the fields, including pre/post_processing,
fknob_creator/deletor, matrices_to_track, reordering_power, and reordering_FAR
are all for the new first statement.

One special case frequently seen from a source line like x += ... is:
    y = f(...) # f is some function. y is a Symbol or GenSym
    result = y
We would replace the two statements into the following instead:
    f!(result, ...)
    y = result
Result must have already been allocated space, and that space is not used by
any other array. Otherwise, we cannot do the replacement.
ASSUMPTION: the only difference between f and f! is that f! has the result 
as its first parameter (Space already allocated), while f must allocate space
for its result before it returns. By changing f to f!, we save memory allocation.
"""
immutable TwoStatementsPattern <: Pattern
    name                :: AbstractString # Name of the pattern
    first_skeleton      :: Tuple  # Skeleton of the first expression
    second_skeleton     :: Tuple  # Skeleton of the second expression
    pre_processing      :: Function
    new_first_skeleton  :: Tuple  # The skeleton for the substitute of the first expression
    new_second_skeleton :: Tuple  # The skeleton for the substitute of the second expression
    post_processing     :: Function
    fknob_creator       :: AbstractString
    fknob_deletor       :: AbstractString
    matrices_to_track   :: Tuple
    reordering_power    :: Int
    reordering_FAR      :: Tuple
end

# Below are the expr_patterns we care about.

# Patterns that are used only for matching (sub-expressions).
@doc """ SpMV(a, A, x) """
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
    (),
    0,
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
    (),
    0,
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
    (),
    0,
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
    (),
    0,
    ()
)

const number_times_matrix_pattern = ExprPattern(
    "number_times_matrix_pattern",
    (:call, GlobalRef(Main, :*), Number, SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE,),
    do_nothing,
    "",
    "",
    (),
    0,
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
    (),
    0,
    ()
)

const vector_minus_number_pattern = ExprPattern(
    "vector_minus_number_pattern",
    (:call, GlobalRef(Main, :-), Vector, Number),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

const vector_add_number_pattern = ExprPattern(
    "vector_add_number_pattern",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    (),
    0,
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
    (),
    0,
    ()
)

# Patterns that do transformation
@doc """ dot(x, y) ==> SparseAccelerator.dot(x, y)"""
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
    (),
    0,
    ()
)

@doc """ norm(x) ==> SparseAccelerator.norm(x)"""
const norm_pattern1 = ExprPattern(
    "norm_pattern1",
    (:call, GlobalRef(Main, :norm), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:norm)),
     :arg2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ sum(x) ==> SparseAccelerator.sum(x)"""
const sum_pattern1 = ExprPattern(
    "sum_pattern1",
    (:call, GlobalRef(Main, :sum), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:sum)),
     :arg2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ mean(x) ==> SparseAccelerator.sum(x)/length(x)"""
const mean_pattern1 = ExprPattern(
    "mean_pattern1",
    (:call, GlobalRef(Main, :mean), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Main, :(/)),
      TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:sum), :arg2),
      TypedExprNode(Function, :call, TopNode(:getfield), :Base, QuoteNode(:arraylen), :arg2)
    ),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ minimum(x) ==> SparseAccelerator.minimum(x)"""
const minimum_pattern1 = ExprPattern(
    "minimum_pattern1",
    (:call, GlobalRef(Main, :minimum), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:minimum)),
     :arg2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ abs!(w, x) ==> SparseAccelerator.abs!(w, x)"""
const abs!_pattern1 = ExprPattern(
    "abs!_pattern1",
    (:call, GlobalRef(Main, :abs!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:abs!)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ exp!(w, x) ==> SparseAccelerator.exp!(w, x)"""
const exp!_pattern1 = ExprPattern(
    "exp!_pattern1",
    (:call, GlobalRef(Main, :exp!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:exp!)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ log1p!(w, x) ==> SparseAccelerator.log1p!(w, x)"""
const log1p!_pattern1 = ExprPattern(
    "log1p!_pattern1",
    (:call, GlobalRef(Main, :log1p!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:log1p!)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ min!(w, x) ==> SparseAccelerator.min!(w, x)"""
const min!_pattern1 = ExprPattern(
    "min!_pattern1",
    (:call, GlobalRef(Main, :min!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:min!)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ copy!(y, x) ==> SparseAccelerator.copy!(y, x)"""
const copy!_pattern1 = ExprPattern(
    "copy!_pattern1",
    (:call, GlobalRef(Main, :copy!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:copy!)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ A * x => SpMV(A, x) """
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
    (),
    0,
    ()
)

@doc """ a * A * x + g => SpMV(a, A, x, 0, x, g) """
const SpMV_pattern2 = ExprPattern(
    "SpMV_pattern2",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, number_times_matrix_vector_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ a * A * x => SpMV(a, A, x) """
const SpMV_pattern3 = ExprPattern(
    "SpMV_pattern3",
    (:call, GlobalRef(Main, :*), Number, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3, :arg4),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ SpMV(a, A, x) + g => SpMV(a, A, x, 0, x, g) """
const SpMV_pattern4 = ExprPattern(
    "SpMV_pattern4",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, SpMV_3_parameters_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ A_mul_B!(y, A, x) = SpMV!(y, A, x) """
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
    (),
    0,
    ()
)

@doc """ z = SpMV(a, A, x, b, y, g), z is x or y => SpMV!(z, a, A, x, b, y, g) """
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
    (),
    0,
    ()
)

@doc """ x = a * A * x => SpMV!(x, a, A, x, 0, x, 0) """
const SpMV!_pattern3 = ExprPattern(
    "SpMV!_pattern3",
    (:(=), Vector, Vector),
    (nothing, nothing, SpMV_3_parameters_pattern),
    LHS_in_RHS,
    (TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg1, :aarg22, :aarg23, :aarg24, 0.0, :arg1, 0.0),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ A_mul_B!(a, A, x, b, y) => SpMV!(y, a, A, x, b, y , 0) """
const SpMV!_pattern4 = ExprPattern(
    "SpMV_pattern4",
    (:call, GlobalRef(Main, :A_mul_B!), Number, SparseMatrixCSC, Vector, Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg6, :arg2, :arg3, :arg4, :arg5, :arg6, 0.0),
    do_nothing,
    "",
    "",
    (),
    0,
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
    (),
    0,
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
    (),
    0,
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
    (),
    0,
    ()
)

const WAXPB!_pattern1 = ExprPattern(
    "WAXPB!_pattern1",
    (:(=), Vector, Vector),
    (nothing, nothing, vector_minus_number_pattern),
    LHS_in_RHS,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPB!)),
     :arg1, 1, :arg1, :naarg23),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

const WAXPB!_pattern2 = ExprPattern(
    "WAXPB!_pattern2",
    (:(=), Vector, Vector),
    (nothing, nothing, vector_add_number_pattern),
    LHS_in_RHS,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPB!)),
     :arg1, 1, :arg1, :aarg23),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ w = x.*y => w = element_wise_multiply(x, y) """
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
    (),
    0,
    ()
)

@doc """ w = x ./ y => w = element_wise_divide(x, y) """
const element_wise_divide_pattern1 = ExprPattern(
    "element_wise_divide_pattern1",
    (:call, GlobalRef(Main, :./), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:element_wise_divide)),
     :arg2, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ Main.trace(A) => SparseAccelerator.trace(A) """
const trace_pattern1 = ExprPattern(
    "trace_pattern1",
    (:call, GlobalRef(Main, :trace), SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:trace)),
     :arg2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
a * A - B => SpAdd(a, A, -1, B)
TODO: make the pattern more general to cover cases like a * A +- b * B. This 
requires a pattern to handle more than one shape (e.g. b = 1 and b != 1 cases,
and + and - cases)
"""
const SpAdd_pattern1 = ExprPattern(
    "SpAdd_pattern1",
    (:call, GlobalRef(Main, :-), SparseMatrixCSC{Float64, Int32}, SparseMatrixCSC{Float64, Int32}),
    (nothing, nothing,           number_times_matrix_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpAdd)),
     :aarg22, :aarg23, -1, :arg3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

expr_patterns = [
    dot_pattern1,
    norm_pattern1,
    sum_pattern1,
    #mean_pattern1, #ISSUE: this pattern requires replacement of arguments in a sub-tree. #TODO: enalble matching and replacing a tree with multiple levels.
    minimum_pattern1, 
    abs!_pattern1, 
    exp!_pattern1, 
    log1p!_pattern1, 
    min!_pattern1,
    copy!_pattern1,
    #WAXPBY_pattern,
    SpMV_pattern1,
    SpMV_pattern2,
    SpMV_pattern3,
    SpMV_pattern4,
    SpMV!_pattern1,
    SpMV!_pattern2,
    SpMV!_pattern3,
    SpMV!_pattern4,
    WAXPBY_pattern1,
    WAXPBY_pattern2,
    WAXPBY!_pattern1,
    WAXPB!_pattern1,
    WAXPB!_pattern2,
    element_wise_multiply_pattern1,
    element_wise_divide_pattern1,
    trace_pattern1,
    SpAdd_pattern1
    #WAXPBY!_pattern,
]

@doc """
    z = SparseAccelerator.SpMV(a, A, x, b, y, r)
    y = z
=>
    SparseAccelerator.SpMV!(y, a, A, x, b, y, r)
    z = y    
"""
const SpMV!_two_statements_pattern1 = TwoStatementsPattern(
    "SpMV!_two_statements_pattern1",
    (:(=), Any, 
           Expr(:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
                 Number, SparseMatrixCSC, Vector, Number, Vector, Number)),
    (:(=), Vector, :f1),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :s1, :f2_2, :f2_3, :f2_4, :f2_5, :f2_6, :f2_7),
    (:(=), :f1, :s1),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    z = SparseAccelerator.WAXPBY(a, x, b, y)
    y = z
=>
    SparseAccelerator.WAXPBY!(y, a, x, b, y)
    z = y    
"""
const WAXPBY!_two_statements_pattern1 = TwoStatementsPattern(
    "WAXPBY!_two_statements_pattern1",
    (:(=), Any,
           Expr(:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY)),
                 Number, Vector, Number, Vector)),
    (:(=), Vector, :f1),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:WAXPBY!)),
     :s1, :f2_2, :f2_3, :f2_4, :f2_5),
    (:(=), :f1, :s1),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

two_statements_patterns = [
    SpMV!_two_statements_pattern1,
    WAXPBY!_two_statements_pattern1
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
Use the real argument in args to replace the symbolic arg (like :e1).
First/second_expr are used only when replacing a two-statements pattern.
"""
function replacement_arg(
    arg         :: Any,    # Symbolic arg to replace
    args        :: Vector, # Real arguments
    result      :: Any,
    symbol_info :: Sym2TypeMap,
    first_expr  :: Any = nothing,
    second_expr :: Any = nothing
)
    if typeof(arg) == Symbol
        arg_string = string(arg)
        if length(arg_string) == 6 && arg_string[1 : 6] == "result" 
            # This refers to the result arg of expr.
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
        elseif length(arg_string) > 1 && (arg_string[1] == 'f' || arg_string[1] == 's')
            # f/s means first/second_expr.
            #   f3 means   first_expr.args[3].
            #   f3_1 means first_expr.args[3].args[1]
            #   s3 means   second_expr.args[3].
            #   s3_1 means second_expr.args[3].args[1]
            # TODO: Get rid of "arg" and "aarg", and replace them all with f or s.
            assert(first_expr  != nothing)
            assert(second_expr != nothing)
            indexes = split(arg_string[2 : end], "_")
            x = parse(Int, indexes[1])
            if arg_string[1] == 'f'
                arg = first_expr.args[x]
            else
                arg = second_expr.args[x]
            end 
            for i in 2 : length(indexes)
                x = parse(Int, indexes[i])
                arg = arg.args[x]
            end         
        end
    end
    arg
end

@doc """
Replace the AST based on the substitute skeleton. Return true if AST is replaced.
When repalcing an ExprPattern, expr is a copy of first_expr, and second_expr is
nothing. When replacing a TwoStatementsPattern, expr is a copy of either
first or second_expr. 
"""
function replace(
    substitute                 :: Tuple,
    expr                       :: Expr,
    symbol_info                :: Sym2TypeMap,
    first_expr                 :: Any,
    second_expr                :: Any,
    expr_is_copy_of_first_expr :: Bool
)
    if length(substitute) == 1 && substitute[1] == :NO_CHANGE
        return false
    end

    expr.head = substitute[1]
    empty!(expr.args)
    for i = 2 : length(substitute)
        arg = replacement_arg(substitute[i], 
                expr_is_copy_of_first_expr ? first_expr.args : second_expr.args,
                nothing, symbol_info, first_expr, second_expr)
        push!(expr.args, arg)
    end
    return true
end

function replace(
    substitute  :: Tuple,
    expr        :: Expr,
    symbol_info :: Sym2TypeMap,
) 
    expr_backup  = copy(expr)
    replace(substitute, expr, symbol_info, expr_backup, nothing, true)
end

@doc """ 
Match the expr's skeleton with the pattern's skeleton. Return true if matched.
First/second_expr are used only when matching a two-statements pattern, where
expr is part or whole of first/second_expr. 
"""
function match_skeletons(
    expr             :: Expr,
    pattern_skeleton :: Tuple,
    symbol_info      :: Sym2TypeMap,    
    first_expr       :: Any = nothing,
    second_expr      :: Any = nothing,          
)
    skeleton = expr_skeleton(expr, symbol_info)
    if length(skeleton) == length(pattern_skeleton)
        if (skeleton[1] == pattern_skeleton[1])
            # So far, we handle only call or assignment
            assert((skeleton[1] == :call) || (skeleton[1] == :(=)))
            for i in 2 : length(skeleton)
                # A skeleton has a head in the first place, then arguments.
                # Thus expr.args[i - 1] corresponds to pattern_skeleton[i]
                real_arg = expr.args[i - 1]      

                typ = typeof(pattern_skeleton[i]) 
                if typ == Expr
                    if pattern_skeleton[i].typ == Function
                        # Literally match the module and function name
                        if skeleton[i] != pattern_skeleton[i]
                            return false
                        else
                            continue
                        end
                    end
                    
                    if typeof(real_arg) != Expr
                        return false
                    else
                        # Do recursive match
                        sub_pattern_skeleton = ntuple(j -> (
                            (j == 1) ? pattern_skeleton[i].head
                                     : pattern_skeleton[i].args[j - 1]),
                            length(pattern_skeleton[i].args) + 1)
                        if !match_skeletons(real_arg, sub_pattern_skeleton, symbol_info, first_expr, second_expr)
                            return false
                        end
                    end
                elseif typ == Symbol
                    # pattern_skeleton[i] is a symbolic argument like :f1 or :s1.
                    arg = replacement_arg(pattern_skeleton[i], expr.args, nothing, symbol_info, first_expr, second_expr)
                    if get_symexpr(arg) != get_symexpr(real_arg)
                        return false
                    end
                elseif typ == GlobalRef
                    # We need literal match in these cases.
                    if skeleton[i] != pattern_skeleton[i]
                        return false
                    end
                elseif !(skeleton[i] <: pattern_skeleton[i])
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
    if match_skeletons(ast, pattern.skeleton, symbol_info)
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
                        if !match_replace(sub_pattern, ast.args[i - 1], call_sites)
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

        dprintln(1, 1, "Matched ", pattern.name, " with ", ast)
    
        ast_changed = replace(pattern.substitute, ast, symbol_info)

        if ast_changed
            dprintln(1, 2, "Replaced with ", ast)
        end

        if pattern.post_processing != do_nothing
            if !pattern.post_processing(ast, call_sites, pattern.fknob_creator,
                                        pattern.fknob_deletor, pattern.matrices_to_track,
                                        pattern.reordering_power, pattern.reordering_FAR)
                # AST has already been changed by replace(). However, post 
                # processing fails. That AST might be wrong. So abort 
                dprintln(1, 2, "Post-processing failed.")
                throw(PostPatternReplacementFailure(pattern))
            end
        end

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
        for pattern in patterns
            success = match_replace(pattern, ast, call_sites)
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

@doc """
Check if the two expression, from adjacent statements, are of the following pattern:
    y = f(...) # f is some function. y is a Symbol or GenSym
    result = y
Result must have already been allocated space.

Replace it into the following instead:
    f!(result, ...)
    y = result
"""
function match_replace_an_two_statements_pattern(
    first_expr  :: Expr,
    second_expr :: Expr,
    call_sites  :: CallSites,
    patterns    :: Vector
)      
    # Check that the result is one of the input parameter of f, in which case
    # it must have been allocated space already.
    # ISSUE: this does not ensure safety, because the result's memory
    # might be shared(pointed to) by another array, and thus we cannot
    # do this replacement.
    # NOTE: In a for loop, Julia'a semantics is to create a new 
    # result array every time result = ... is encountered.
    # TODO: after points-to and memory lifetime analysis work, check that
    # the memory is not only allocated, but also no other array uses it.  
    #if !is_an_arg(expr.args[1], prev_expr.args[2])
    #    return false
    #end
    symbol_info = call_sites.symbol_info    
    for pattern in patterns
        if match_skeletons(first_expr,  pattern.first_skeleton,  symbol_info, first_expr, second_expr) &&
           match_skeletons(second_expr, pattern.second_skeleton, symbol_info, first_expr, second_expr)
            dprintln(1, 1, "Matched ", pattern.name, " with: ")
            dprintln(1, 2, first_expr)
            dprintln(1, 2, second_expr)

            if pattern.pre_processing != do_nothing
                if !pattern.pre_processing(first_expr, call_sites, pattern.fknob_creator,
                                           pattern.fknob_deletor)
                    return false
                end
            end

            first_expr_backup  = copy(first_expr)
            second_expr_backup = copy(second_expr)
            ast_changed        = replace(pattern.new_first_skeleton, first_expr, symbol_info, first_expr_backup, second_expr_backup, true)
            ast_changed       |= replace(pattern.new_second_skeleton, second_expr, symbol_info, first_expr_backup, second_expr_backup, false)

            if pattern.post_processing != do_nothing
                if !pattern.post_processing(first_expr, call_sites, pattern.fknob_creator,
                                            pattern.fknob_deletor, pattern.matrices_to_track,
                                            pattern.reordering_power, pattern.reordering_FAR)
                    # AST has already been changed by replace(). However, post 
                    # processing fails. That AST might be wrong. So abort 
                    dprintln(1, 2, "Post-processing failed.")
                    throw(PostPatternReplacementFailure(pattern))
                end
            end

            if ast_changed
                dprintln(1, 2, "Replaced with:")
                dprintln(1, 2, first_expr)
                dprintln(1, 2, second_expr)
            end

            return true
        end
    end
    return false
end

@doc """ 
Pattern match and replace the code that is functionally equivalent to SpMV, SpMV!,
dot, WAXPBY, WAXPBY!, etc. with calls to the corresponding SPMP library functions.
"""
function replace_calls(
    func_ast    :: Expr,
    symbol_info :: Sym2TypeMap, 
    cfg         :: CFG
)
    dprintln(1, 0, "\nReplacing sparse matrix function calls:")

    call_sites  = CallSites(Set{CallSite}(), FunctionRegion(func_ast), symbol_info, 
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
                # Try to merge this and the previous expression
                match_replace_an_two_statements_pattern(prev_expr, expr, call_sites, two_statements_patterns)
            end
            
            prev_expr = expr
        end
    end
end
