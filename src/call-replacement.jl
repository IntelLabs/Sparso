#=
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

TypedExprNode(x...) = LivenessAnalysis.TypedExpr(x...)

@doc """
The CallSites' extra field for call replacement.
"""
type CallReplacementExtra
    # Environment before the AST walking
    matrix_properties        :: Symexpr2PropertiesMap

    # Temporary variables created to represent the result of an expression.
    expression_temporaries   :: Set{Sym}
    
    # Scratch variables.
    live_in_before_prev_expr :: Set{Sym}
    live_in_before_expr      :: Set{Sym}
    bb                       :: Any              # BasicBlock
    stmt_idx                 :: StatementIndex

    # Some patterns (those contain :t1!2, etc.) may hoist some subtrees of the
    # current statement before it. That splits one statement into more than one.
    # Such patterns should be matched at the last, because otherwise, other 
    # patterns may not be able to match what they should: they cannot find the
    # subtrees to match, which are no longer in the same statement.
    # We call such patterns splitting patterns, and the other non-splitting.
    # We can top-down match non-splitting patterns once, and bottom-up match
    # all patterns (including splitting and non-splitting) once.
    non_splitting_patterns  :: Vector{Pattern}
        
    CallReplacementExtra(_matrix_properties) = new(_matrix_properties, 
            Set{Sym}(), Set{Sym}(), Set{Sym}(), nothing, 0, Vector{Pattern}())
end

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


abstract AbstractPatternAction

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
(the pattern is for matching only); :a1 means the Expr's args[1], :a1_2 means
the Expr's args[1].args[2], :n1_2 means the negative :a1_2, etc.

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
    pre_processing    :: Union{Function, AbstractPatternAction}
    substitute        :: Tuple
    post_processing   :: Union{Function, AbstractPatternAction}
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

@doc """
Argument Description. Used for describing an argument in a pattern.
type_or_symbol: a type of, or a symbol for, an argument.
properties:     a combination of properties. 
relations:      relation with other argument/variable.
 
For example,
    AD(SparseMatrixCSC, SA_HAS_FREE_MEMORY | SA_CONST_VALUED, [(SA_TRANSPOSE_OF, :a1])
means: the current argument has SparseMatrixCSC type, it has free memory
available at the program point the argument is executed, it is 
is constant valued (in the current region), and it is a transposition of
args[1] in the same expression.

For another example,
    AD(:f2_3, SA_SYMM_STRUCTURED) 
means that in a TwoStatementPattern, an argument that is the same as the 
first_expr.args[2].args[3], and it is symmetric in structure.

A symbol can be simple or more complicated like :a1_2, :f1_2, :s1_2, :n1, etc.
For example:
    :a1 means the current Expr's args[1], 
    :a1_2 means the current Expr's args[1].args[2], 
    :n1_2 means the negative :a1_2
    :f1_2 means (in a TwoStatementPattern) first_expr.args[1].arg[2]
    :s1_2 means (in a TwoStatementPattern) second_expr.args[1].arg[2]

When you use AD, you need to describe at least properties or relations. 
Otherwise, you can directly describe it in a type or a symbol, without using
AD.
"""
typealias Property Int                        # e.g. SA_CONST_STRUCTURED | SA_SYMM_STRUCTURED
typealias Relation Tuple{Property, Property}  # e.g. (SA_UPPER_OF, SA_CONST_STRUCTURED | SA_SYMM_STRUCTURED), 
                                              # which means the structure of the first matrix is
                                              # upper of that of the second matrix, and the second
                                              # matrix is constant and symmetric in structure 
type AD
    type_or_symbol :: Any
    properties     :: Vector{Union{Property, Relation}}
    
    AD(_type_or_symbol) = new(_type_or_symbol, Vector{Union{Property, Relation}}())
    function AD(_type_or_symbol, _properties...) 
        ad = new(_type_or_symbol, Vector{Union{Property, Relation}}())
        for property in _properties
            push!(ad.properties, property)
        end
        ad
    end
end

# Below are the expr_patterns we care about.

# Patterns that are used only for matching (sub-expressions).
@doc """ SpMV(a, A, x) """
const SpMV_3_parameters_pattern = ExprPattern(
    "SpMV_3_parameters_pattern",
    (:call, GlobalRef(Sparso, :SpMV),
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
    (:call, GlobalRef(Sparso, :SpMV),
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
    (:call, GlobalRef(Sparso, :SpMV),
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

@doc """ x - a """
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

@doc """ x + a """
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

@doc """ Sparso.WAXPBY(a, x, b, y) """
const WAXPBY_4_parameters_pattern = ExprPattern(
    "WAXPBY_4_parameters_pattern",
    (:call, GlobalRef(Sparso, :WAXPBY),
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
@doc """ dot(x, y) ==> Sparso.dot(x, y)"""
const dot_pattern1 = ExprPattern(
    "dot_pattern1",
    (:call, GlobalRef(Main, :dot), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :dot), :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ norm(x)^2 ==> Sparso.dot(x, x)"""
const dot_pattern2 = ExprPattern(
    "dot_pattern2",
    (:call, GlobalRef(Main, :^), Expr(:call, GlobalRef(Main, :norm), Vector), 2),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :dot),
     :a2_2, :a2_2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ Sparso.norm(x)^2 ==> Sparso.dot(x, x)"""
const dot_pattern3 = ExprPattern(
    "dot_pattern3",
    (:call, GlobalRef(Main, :^), Expr(:call, GlobalRef(Sparso, :norm), Vector), 2),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :dot),
     :a2_2, :a2_2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ norm(x) ==> Sparso.norm(x)"""
const norm_pattern1 = ExprPattern(
    "norm_pattern1",
    (:call, GlobalRef(Main, :norm), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :norm), :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ sum(x) ==> Sparso.sum(x)"""
const sum_pattern1 = ExprPattern(
    "sum_pattern1",
    (:call, GlobalRef(Main, :sum), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :sum),
     :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ mean(x) ==> Sparso.sum(x)/length(x)"""
const mean_pattern1 = ExprPattern(
    "mean_pattern1",
    (:call, GlobalRef(Main, :mean), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Main, :(/)),
      (:call, GlobalRef(Sparso, :sum), :a2),
      (:call, GlobalRef(Base, :arraylen), :a2)
    ),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ minimum(x) ==> Sparso.minimum(x)"""
const minimum_pattern1 = ExprPattern(
    "minimum_pattern1",
    (:call, GlobalRef(Main, :minimum), Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :minimum),
     :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ abs!(w, x) ==> Sparso.abs!(w, x)"""
const abs!_pattern1 = ExprPattern(
    "abs!_pattern1",
    (:call, GlobalRef(Main, :abs!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :abs!),
     :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    abs(x) ==> Sparso.abs!(temp, x), where temp is a temporary to be
    generated by the compiler.
"""
const abs!_pattern2 = ExprPattern(
    "abs!_pattern2",
    (:call, GlobalRef(Main, :abs), Vector{Float64}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :abs!),
     :t2, :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ exp!(w, x) ==> Sparso.exp!(w, x)"""
const exp!_pattern1 = ExprPattern(
    "exp!_pattern1",
    (:call, GlobalRef(Main, :exp!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :exp!),
     :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    exp(x) ==> Sparso.exp!(temp, x), where temp is a temporary to be
    generated by the compiler.
"""
const exp!_pattern2 = ExprPattern(
    "exp!_pattern2",
    (:call, GlobalRef(Main, :exp), Vector{Float64}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :exp!),
     :t2, :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    log(1 + x) ==> Sparso.log1p!(temp, x), where temp is a temporary to be
    generated by the compiler.
"""
const log1p!_pattern1 = ExprPattern(
    "log1p!_pattern1",
    (:call, GlobalRef(Main, :log), Expr(:call, GlobalRef(Main, :+), 1, Vector{Float64})),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :log1p!),
     :t2_3, :a2_3), # If use :t2 instead of :t2_3, since arg2 is an Expr and thus will be hoisted ahead of the statement and replaced with a symbol, the next :a2_3 will be meaningless.
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    log(Sparso.WAXPB!(temp, 1, x, 1) ==> Sparso.log1p!(temp, x), 
    where temp is a temporary already generated by the compiler. We can use
    temp as the result of log because a temp is only used in the tree and only once,
    so that we are sure it will not be used outside the tree, nor used by any other 
    nodes in the tree.
"""
const log1p!_pattern2 = ExprPattern(
    "log1p!_pattern2",
    (:call, GlobalRef(Main, :log), 
      Expr(:call, GlobalRef(Sparso, :WAXPB!), Vector{Float64}, 1, Vector{Float64}, 1)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :log1p!),
     :t2_2, :a2_4), # If use :t2 instead of :t2_2, since arg2 is an Expr and thus will be hoisted ahead of the statement and replaced with a symbol, the next :a2_4 will be meaningless.
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    min(x, a) ==> Sparso.min!(temp, x, a), where temp is a temporary to be
    generated by the compiler.
"""
const min!_pattern1 = ExprPattern(
    "min!_pattern1",
    (:call, GlobalRef(Main, :min), Vector{Float64}, Number),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :min!),
     :t2, :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ copy!(y, x) ==> Sparso.copy!(y, x) """
const copy!_pattern1 = ExprPattern(
    "copy!_pattern1",
    (:call, GlobalRef(Main, :copy!), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :copy!), :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    copy(x) ==> Sparso.copy!(temp, x), where temp is a temporary to be
    generated by the compiler.
"""
const copy!_pattern2 = ExprPattern(
    "copy!_pattern2",
    (:call, GlobalRef(Main, :copy), Vector{Float64}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :copy!), :t0!2, :a2),
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
    (:call, GlobalRef(Sparso, :SpMV),
     :a2, :a3),
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
    (:call, GlobalRef(Sparso, :SpMV),
     :a2_2, :a2_3, :a2_4, 0, :a2_4, :a3),
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
    (:call, GlobalRef(Sparso, :SpMV),
     :a2, :a3, :a4),
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
    (:call, GlobalRef(Sparso, :SpMV),
     :a2_2, :a2_3, :a2_4, 0, :a2_4, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """
    -(A*x)/a + b*y => SpMV(-1/a, A, x, b, y)
TODO: check why AST walker matches subtrees first, and make it match a parent
tree first. That would make his pattern work.
"""
const SpMV_pattern5 = ExprPattern(
    "SpMV_pattern5",
    (:call, GlobalRef(Main, :+), 
      Expr(:call, GlobalRef(Main, :/), 
        Expr(:call, GlobalRef(Main, :-),      
              Expr(:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector)),
        Number),     
      Expr(:call, GlobalRef(Main, :*), Number, Vector)), 
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV),
     Expr(:call, GlobalRef(Main, :/), -1, :a2_3), :a2_2_2_2, :a2_2_2_3, :a3_2, :a3_3),
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
    (:call, GlobalRef(Sparso, :SpMV!),
     :a2, :a3, :a4),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ z = SpMV(a, A, x, b, y, g) => SpMV!(z, a, A, x, b, y, g) """
const SpMV!_pattern2 = ExprPattern(
    "SpMV!_pattern2",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), Vector),
    (nothing, nothing, SpMV_6_parameters_pattern),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :a1, :a2_2, :a2_3, :a2_4, :a2_5, :a2_6, :a2_7),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ y = a * A * x => SpMV!(y, a, A, x, 0, x, 0) """
const SpMV!_pattern3 = ExprPattern(
    "SpMV!_pattern3",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), Vector),
    (nothing, nothing, SpMV_3_parameters_pattern),
    do_nothing,
    (GlobalRef(Sparso, :SpMV!),
     :a1, :a2_2, :a2_3, :a2_4, 0.0, :a1, 0.0),
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
    (:call, GlobalRef(Sparso, :SpMV!),
     :a6, :a2, :a3, :a4, :a5, :a6, 0.0),
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
    (:call, GlobalRef(Sparso, :WAXPBY),
     1, :a2, :a3_2, :a3_3),
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
    (:call, GlobalRef(Sparso, :WAXPBY),
     1, :a2, :n3_2, :a3_3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ z = Sparso.WAXPBY(a, x, b, y) => Sparso.WAXPBY!(z, a, x, b, y) """
const WAXPBY!_pattern1 = ExprPattern(
    "WAXPBY!_pattern1",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), 
           Expr(:call, GlobalRef(Sparso, :WAXPBY), Number, Vector, Number, Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPBY!),
     :a1, :a2_2, :a2_3, :a2_4, :a2_5),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ 
x - y => Sparso.WAXPBY!(temp, 1, x, -1, y) , where temp is a temporary
to be generated by the compiler.
"""
const WAXPBY!_pattern2 = ExprPattern(
    "WAXPBY!_pattern2",
    (:call, GlobalRef(Main, :-), Vector{Float64}, Vector{Float64}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPBY!),
     :t2!3, 1, :a2, -1, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ 
 -x => Sparso.WAXPBY!(temp, -1, x, 0, x) , where temp is a temporary
to be generated by the compiler.
"""
const WAXPBY!_pattern3 = ExprPattern(
    "WAXPBY!_pattern3",
    (:call, GlobalRef(Main, :-), Vector{Float64}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPBY!),
     :t2, -1, :a2, 0, :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ w = x - a => Sparso.WAXPB!(w, 1, x, -a)"""
const WAXPB!_pattern1 = ExprPattern(
    "WAXPB!_pattern1",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), Vector),
    (nothing, nothing, vector_minus_number_pattern),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPB!),
     :a1, 1, :a2_2, :n2_3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ w = x + a => Sparso.WAXPB!(w, 1, x, a)"""
const WAXPB!_pattern2 = ExprPattern(
    "WAXPB!_pattern2",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), Vector),
    (nothing, nothing, vector_add_number_pattern),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPB!),
     :a1, 1, :a2_2, :a2_3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ 
a + x => Sparso.WAXPB!(temp, 1, x, a), where temp is a temporary
to be genreated by the compiler.
"""
const WAXPB!_pattern3 = ExprPattern(
    "WAXPB!_pattern3",
    (:call, GlobalRef(Main, :+), Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPB!),
     :t3, 1, :a3, :a2),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ x.*y => w = element_wise_multiply(x, y) """
const element_wise_multiply_pattern1 = ExprPattern(
    "element_wise_multiply_pattern1",
    (:call, GlobalRef(Main, :.*), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :element_wise_multiply),
     :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ w = x .* y => Sparso.element_wise_multiply!(w, x, y)"""
const element_wise_multiply!_pattern1 = ExprPattern(
    "element_wise_multiply!_pattern1",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), 
        Expr(:call, GlobalRef(Sparso, :element_wise_multiply), Vector, Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :element_wise_multiply!),
     :a1, :a2_2, :a2_3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ x ./ y => element_wise_divide(x, y) """
const element_wise_divide_pattern1 = ExprPattern(
    "element_wise_divide_pattern1",
    (:call, GlobalRef(Main, :./), Vector, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :element_wise_divide),
     :a2, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ w = x ./ y => Sparso.element_wise_divide!(w, x, y)"""
const element_wise_divide!_pattern1 = ExprPattern(
    "element_wise_divide!_pattern1",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), 
        Expr(:call, GlobalRef(Sparso, :element_wise_divide), Vector, Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :element_wise_divide!),
     :a1, :a2_2, :a2_3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

@doc """ Main.trace(A) => Sparso.trace(A) """
const trace_pattern1 = ExprPattern(
    "trace_pattern1",
    (:call, GlobalRef(Main, :trace), SparseMatrixCSC),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :trace),
     :a2),
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
    (:call, GlobalRef(Sparso, :SpAdd),
     :a2_2, :a2_3, -1, :a3),
    do_nothing,
    "",
    "",
    (),
    0,
    ()
)

expr_patterns = [
    dot_pattern1,
    dot_pattern2,
    dot_pattern3,
    norm_pattern1,
    sum_pattern1,
    #mean_pattern1, #ISSUE: this pattern requires replacement of arguments in a sub-tree. #TODO: enalble matching and replacing a tree with multiple levels.
    minimum_pattern1, 
    abs!_pattern1,
    abs!_pattern2, 
    exp!_pattern1,
    exp!_pattern2,
    log1p!_pattern1, 
    log1p!_pattern2, 
    min!_pattern1,
    copy!_pattern1,
    copy!_pattern2,
    #WAXPBY_pattern,
    SpMV_pattern5,
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
    WAXPBY!_pattern2,
    WAXPBY!_pattern3,
    WAXPB!_pattern1,
    WAXPB!_pattern2,
    WAXPB!_pattern3,
    element_wise_multiply_pattern1,
    element_wise_multiply!_pattern1,
    element_wise_divide_pattern1,
    element_wise_divide!_pattern1,
    trace_pattern1,
    SpAdd_pattern1
    #WAXPBY!_pattern,
]

@doc """
    z = Sparso.SpMV(a, A, x, b, y, r)
    y = z
=>
    Sparso.SpMV!(y, a, A, x, b, y, r)
    z = y    
"""
const SpMV!_two_statements_pattern1 = TwoStatementsPattern(
    "SpMV!_two_statements_pattern1",
    (:(=), Any, 
           Expr(:call, GlobalRef(Sparso, :SpMV),
                 Number, SparseMatrixCSC, Vector, Number, Vector, Number)),
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), :f1),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
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
    z = Sparso.WAXPBY(a, x, b, y)
    y = z
=>
    Sparso.WAXPBY!(y, a, x, b, y)
    z = y    
"""
const WAXPBY!_two_statements_pattern1 = TwoStatementsPattern(
    "WAXPBY!_two_statements_pattern1",
    (:(=), Any,
           Expr(:call, GlobalRef(Sparso, :WAXPBY),
                 Number, Vector, Number, Vector)),
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY), :f1),
    do_nothing,
    (:call, GlobalRef(Sparso, :WAXPBY!),
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
Use the real argument in args to replace arg, the symbolic argument. Return
the real argument.
First/second_expr are used only when replacing a two-statements pattern.
"""
function replacement_arg(
    arg                 :: Any,    # Symbolic arg to replace
    real_args           :: Vector, # Real arguments
    result              :: Any,
    call_sites          :: CallSites,
    first_expr          :: Any = nothing,
    second_expr         :: Any = nothing,
    trace_pattern_match :: Bool = false              
    
)
    if trace_pattern_match
        println("\tReplacing arg: ", arg, " real_args: ", real_args)
        flush(STDOUT)
    end
    
    symbol_info = call_sites.symbol_info
    if typeof(arg) == Symbol
        arg_string = string(arg)
        if arg == :r 
            # This refers to the result arg of expr.
            assert(result != nothing)
            arg = result
        elseif length(arg_string) > 1 && (arg_string[1] == 't')
            arg = get_temporary(arg_string, real_args, result, call_sites)
        elseif length(arg_string) > 1 && (arg_string[1] == 'f' || arg_string[1] == 's' ||
            arg_string[1] == 'n' || arg_string[1] == 'a')
            # a means "argument"; n means "negative of the argument".
            # f/s means a first/second_expr's argument.
            #   f3 means   first_expr.args[3].
            #   f3_1 means first_expr.args[3].args[1]
            #   s3 means   second_expr.args[3].
            #   s3_1 means second_expr.args[3].args[1]
            assert((arg_string[1] != 'f' && arg_string[1] != 's')|| first_expr  != nothing)
            assert((arg_string[1] != 'f' && arg_string[1] != 's')|| second_expr != nothing)
            indexes = split(arg_string[2 : end], "_")
            x = parse(Int, indexes[1])
            if arg_string[1] == 'f'
                arg = first_expr.args[x]
            elseif arg_string[1] == 's'
                arg = second_expr.args[x]
            else
                arg = real_args[x]
            end 
            
            for i in 2 : length(indexes)
                x = parse(Int, indexes[i])
                arg = arg.args[x]
            end         

            if arg_string[1] == 'n'
                arg = TypedExprNode(type_of_ast_node(arg, symbol_info), :call, :(-), arg)
            end
        end
    elseif typeof(arg) == Expr
        new_expr = copy(arg)
        empty!(new_expr.args)
        for a in arg.args
            new_a = replacement_arg(a, real_args, result, call_sites, first_expr, second_expr, trace_pattern_match) 
            push!(new_expr.args, new_a)
        end
        arg = new_expr
    end
    arg
end

@doc """
Get a temporary as specified by a symbolic argument in a pattern.
Arg_string is the string of the symbolic argument, starting with 't'.
"t1_2!3!4" means the argument can be replaced by any of real_args[1].args[2], 
real_args[3], and real_args[4], as long as it is a temporary in an expression.
If none of them is a temporary, a temporary has to be created and allocated
memory, and the corresponding action is recorded.
"t0!1_2!3!4" is special: when a number right after 't' is 0, it means that a
temporary has to be created and allocated, and any of the arguments showing
there can be used to extract the size (and type) info.
"""
function get_temporary(
    arg_string :: AbstractString,
    real_args  :: Vector,
    result              :: Any,    
    call_sites :: CallSites
)
    assert(length(arg_string) > 1 && (arg_string[1] == 't'))
    # So far, :t.. is used only in patterns in the call replacement phase.
    assert(typeof(call_sites.extra) == CallReplacementExtra)

    symbol_info     = call_sites.symbol_info
    possibilities   = split(arg_string[2 : end], "!")
    reuse_temp      = true
    example         = nothing # A symbol based on which to create a new temporary, if necessary

    # An argument, which is a sub expression, might have to be hoisted before
    # the current statement, in order to create a temporary. This hoist canidate
    # is recorded in terms of args[x]. 
    hoist_candidate_args = nothing
    hoist_candidate_x    = 0

    for i in 1 : length(possibilities)
        indexes = split(possibilities[i], "_")
        x       = parse(Int, indexes[1])
        if x == 0
            # "t0" means a temporary has to be created and allocated
            reuse_temp = false
            assert(length(indexes) == 1)
            continue
        end
        args = real_args
        arg  = args[x]
        for j in 2 : length(indexes)
            x = parse(Int, indexes[j])
            args = arg.args
            arg  = args[x]
        end         
        
        # Now the specified argument is found. See if it can reused.
        # So far, "t1_2!3!4" etc. is used only for expressing a vector, and only
        # for Float64 element type
        assert(type_of_ast_node(arg, symbol_info) <: Vector{Float64})
        if in(arg, call_sites.extra.expression_temporaries)
            if reuse_temp
                return arg
            else
                example = get_symexpr(arg)
                break
            end
        else
            s = get_symexpr(arg)
            if typeof(s) <: Sym
                example = s
                break
            else
                assert(typeof(s) == Expr)
                if hoist_candidate_args == nothing
                    hoist_candidate_args = args
                    hoist_candidate_x    = x
                end
            end
        end
    end

    # Create a temporary
    #temp = new_symbol("temp")
    len  = length(call_sites.extra.expression_temporaries)
    temp = Symbol(string("__", "temp", string(len), "__"))
    push!(call_sites.extra.expression_temporaries, temp)

    # Make the new temporary memory allocated.
    lambda   = call_sites.lambda
    bb       = call_sites.extra.bb
    stmt_idx = call_sites.extra.stmt_idx
    
    # TODO: uncomment this statement. And remove the hack
    #   action   = InsertBeforeOrAfterStatement(Vector{Statement}(), true, bb, stmt_idx)
    #   push!(call_sites.actions, action)
    
    # ISSUE: what if the example is NOT live into the loop? How can we make a temporary
    # based on it before the loop?
    # A possible solution: create a dynamic global count so that it creates the
    # temporary right before the statement, but only once: only when the dynamic count is
    # less than the temporary's static index number.
    if example != nothing
        typ               = type_of_ast_node(example, symbol_info)
        symbol_info[temp] = typ 
        addLocalVariable(temp, typ, 2 | 16, lambda) # 2: assigned 16: statically assigned only once
        stmt = Statement(0, Expr(:(=), temp, 
                                  TypedExprNode(typ,:call, :Array, :Cdouble, 
                                        Expr(:call, TopNode(:arraylen), 
                                              SymbolNode(example, typ)))))
                                              
        if typeof(call_sites.region) == FunctionRegion
            action = InsertBeforeOrAfterStatement(Vector{Statement}(), true, bb, stmt_idx)
        else
            # HACK!!! We cannot really insert the generating of a temporary before a loop, 
            # without knowing alias information. This is only for now.
            # TODO: remove the hack.
            assert(typeof(call_sites.region) == LoopRegion)
            action = InsertBeforeLoopHead(Vector{Statement}(), call_sites.region.loop, true)
        end
    else
        # Hoist the candidate, which must be an expression, to be before the
        # statement in order to create a temporary
        # Note: In future, with alias info, we might be able to hoist the candiate
        # as far as possible, even to be before the loop: that becomes loop
        # invariant hoisting.
        assert(hoist_candidate_args != nothing)
        candidate         = hoist_candidate_args[hoist_candidate_x]
        typ               = type_of_ast_node(candidate, symbol_info)
        symbol_info[temp] = typ
        addLocalVariable(temp, typ, 2 | 16, lambda) # 2: assigned 16: statically assigned only once
        
        stmt                                    = Statement(0, TypedExprNode(typ, :(=), temp, candidate))
        hoist_candidate_args[hoist_candidate_x] = SymbolNode(temp, typ)

        action = InsertBeforeOrAfterStatement(Vector{Statement}(), true, bb, stmt_idx)
    end
    push!(call_sites.actions, action)
    push!(action.new_stmts, stmt)
    return temp
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
    call_sites                 :: CallSites,
    first_expr                 :: Any,
    second_expr                :: Any,
    expr_is_copy_of_first_expr :: Bool,
    trace_pattern_match        :: Bool = false              

)
    if length(substitute) == 1 && substitute[1] == :NO_CHANGE
        return false
    end

    expr.head = substitute[1]
    empty!(expr.args)
    for i = 2 : length(substitute)
        arg = replacement_arg(substitute[i],
                expr_is_copy_of_first_expr ? first_expr.args : second_expr.args,
                expr, call_sites, first_expr, second_expr, trace_pattern_match)
        push!(expr.args, arg)
    end
    return true
end

function replace(
    substitute          :: Tuple,
    expr                :: Expr,
    call_sites          :: CallSites,
    trace_pattern_match :: Bool = false              
) 
    expr_backup  = copy(expr)
    replace(substitute, expr, call_sites, expr_backup, nothing, true, trace_pattern_match)
end

@doc """
Check if the argument has the property. This currently happens when in a 
pattern, a skeleton specifies some property to check. 
Live_in_before_expr is the live in set before the first statement in the pattern.
"""
function check_property(
    arg                 :: Any,
    property            :: Property,
    call_sites          :: CallSites,
    live_in_before_expr :: Set{Sym},
)
    # So far, property may include the following cases, if any.
    assert((property & SA_CONST_VALUED != 0) || (property & SA_CONST_STRUCTURED != 0) ||
           (property & SA_SYMM_VALUED != 0) || (property & SA_SYMM_STRUCTURED != 0) ||
           (property & SA_STRUCTURE_ONLY != 0) || (property & SA_HAS_FREE_MEMORY != 0))

    matrix_properties = call_sites.extra.matrix_properties
    if typeof(arg) <: Sym || typeof(arg) == SymbolNode
        arg = get_symexpr(arg)
        if haskey(matrix_properties, arg)
            properties = matrix_properties[arg]
            if property & SA_HAS_FREE_MEMORY != 0
                # ISSUE: this is a hack: even if the arg is live into the pattern,
                # it is still possible that the arg is pointed to by another variable,
                # and thus is not free. Also, even if it does have a free memory, we
                # need to ensure that GC has not collected it (Thus compiler needs to 
                # pin it in some way)
                if !properties.has_dedicated_memory && !in(get_symexpr(arg), live_in_before_expr)
                    return false
                end
            end
            
            if (property & SA_CONST_VALUED != 0) && !properties.constant_valued
                return false
            end

            if (property & SA_CONST_STRUCTURED != 0) && !properties.constant_structured
                return false
            end
                
            if (property & SA_SYMM_VALUED != 0) && !properties.is_symmetric
                return false
            end

            if (property & SA_SYMM_STRUCTURED != 0) && !properties.is_structure_symmetric
                return false
            end

            if (property & SA_STRUCTURE_ONLY != 0) && !properties.is_structure_only
                return false
            end
            return true
        end
    end
    return false
end

@doc """
Check if the argument has the relation. This currently happens when in a 
pattern, a skeleton specifies some relation to check. 
Live_in_before_expr is the live in set before the first statement in the pattern.
"""
function check_relation(
    arg                 :: Any,
    relation            :: Relation,
    call_sites          :: CallSites,
    live_in_before_expr :: Set{Sym},
)
    property1         = relation[1]
    assert((property1 == SA_LOWER_OF) || (property1 == SA_UPPER_OF)) # So far, we handle these two cases
    property2         = relation[2]
    matrix_properties = call_sites.extra.matrix_properties
    if typeof(arg) <: Sym || typeof(arg) == SymbolNode
        arg = get_symexpr(arg)
        if haskey(matrix_properties, arg)
            properties = matrix_properties[arg]
            if (property1 & SA_LOWER_OF != 0)  && (properties.lower_of != nothing) && 
                check_property(properties.lower_of, property2, call_sites, live_in_before_expr)
                return true
            end
            if (property1 & SA_UPPER_OF != 0)  && (properties.upper_of != nothing) && 
                check_property(properties.upper_of, property2, call_sites, live_in_before_expr)
                return true
            end
            return false
        end
    end
    return false
end

@doc """
Check if the argument has the properties. This currently happens when in a 
pattern, a skeleton specifies some properties to check. 
Live_in_before_expr is the live in set before the first statement in the pattern.
"""
function check_properties(
    arg                 :: Any,
    properties          :: Vector{Union{Property, Relation}},
    call_sites          :: CallSites,
    live_in_before_expr :: Set{Sym},
)
    for property in properties
        if typeof(property) == Property
            if !check_property(arg, property, call_sites, live_in_before_expr)
                return false
            end
        else 
            assert(typeof(property) == Relation)
            if !check_relation(arg, property, call_sites, live_in_before_expr)
                return false
            end
        end
    end
    return true
end

@doc """ 
Match the expr's skeleton with the pattern's skeleton. Return true if matched.
Live_in_before_expr is the live in set before the first statement of the pattern.
First/second_expr are used only when matching a two-statements pattern, where
expr is part or whole of first/second_expr. 
"""
function match_skeletons(
    expr                :: Expr,
    pattern_skeleton    :: Tuple,
    call_sites          :: CallSites,
    live_in_before_expr :: Set{Sym},
    first_expr          :: Any = nothing,
    second_expr         :: Any = nothing,
    trace_pattern_match :: Bool = false          
)
    symbol_info = call_sites.symbol_info
    skeleton    = expr_skeleton(expr, symbol_info)
    if trace_pattern_match
        println("\t expr: ", expr)
        println("\t skeleton: ", skeleton)
    end

    if length(skeleton) == length(pattern_skeleton)
        if (skeleton[1] == pattern_skeleton[1])
            # So far, we handle only call or assignment
            assert((skeleton[1] == :call) || (skeleton[1] == :(=)))
            for i in 2 : length(skeleton)
                # A skeleton has a head in the first place, then arguments.
                # Thus expr.args[i - 1] corresponds to pattern_skeleton[i]
                real_arg = expr.args[i - 1]      

                typ = typeof(pattern_skeleton[i])
                
                if trace_pattern_match
                    println("\t skeleton[", i, "]=", skeleton[i])
                    println("\t pattern_skeleton[", i, "]=", pattern_skeleton[i], "::", typ)
                end

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
                        if !match_skeletons(real_arg, sub_pattern_skeleton, call_sites, live_in_before_expr, first_expr, second_expr, trace_pattern_match)
                            return false
                        end
                    end
                elseif typ == Symbol
                    # pattern_skeleton[i] is a symbolic argument like :f1 or :s1.
                    arg = replacement_arg(pattern_skeleton[i], expr.args, nothing, call_sites, first_expr, second_expr, trace_pattern_match)
                    if get_symexpr(arg) != get_symexpr(real_arg)
                        return false
                    end
                elseif typ == AD
                    type_or_symbol = pattern_skeleton[i].type_or_symbol
                    properties     = pattern_skeleton[i].properties
                    #TODO: handle pattern_skeleton[i].relations as well
                    if typeof(type_or_symbol) == Symbol
                        arg = replacement_arg(type_or_symbol, expr.args, expr, call_sites, first_expr, second_expr, trace_pattern_match)
                        if get_symexpr(arg) != get_symexpr(real_arg)
                            return false
                        end
                    else # A type
                        if !(skeleton[i] <: type_or_symbol)
                            return false
                        end
                        arg = real_arg
                    end
                    if !check_properties(arg, properties, call_sites, live_in_before_expr)
                        return false
                    end
                elseif typ <: Number
                    if real_arg != pattern_skeleton[i]
                        return false
                    end
                elseif typ == GlobalRef || typ == TopNode
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

    # For debugging purpose only
    trace_pattern_match = false
    if false #pattern == CS_fwdTriSolve_backslash_pattern && ast.head == :(=) && length(ast.args) ==2 && 
        #typeof(ast.args[2]) == Expr &&
        #ast.args[2].head == :call && ast.args[2].args[1] == GlobalRef(Sparso, :\)
        trace_pattern_match = true
        println("... Matching ", pattern.name, " with ", ast)
    end
        
    live_in          = call_sites.extra.live_in_before_expr
    symbol_info      = call_sites.symbol_info
    match_filter_val = 0 # unknown
    match_filter     = nothing
    if isdefined(call_sites.extra, :pattern_match_filter)
        match_filter = call_sites.extra.pattern_match_filter
        key = tuple(ast, pattern)
        if haskey(match_filter, key)
            match_filter_val = match_filter[key]
            assert(match_filter_val != 0)
        end
    end

    if match_filter_val == 0
        matched = match_skeletons(ast, pattern.skeleton, call_sites, live_in, nothing, nothing, trace_pattern_match)
        if match_filter != nothing
            match_filter[key] = matched ? 1 : -1;
        end
    else
        matched = match_filter_val > 0 ? true : false 
    end

    if matched
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

        #dprintln(1, 1, "Matched ", pattern.name, " with ", ast)
    
        ast_changed = replace(pattern.substitute, ast, call_sites, trace_pattern_match)

        if ast_changed
            dprintln(1, 2, "Replaced with ", ast)
        end

        if pattern.post_processing != do_nothing
            if isa(pattern.post_processing, Function)
                if !pattern.post_processing(ast, call_sites, pattern.fknob_creator,
                                        pattern.fknob_deletor, pattern.matrices_to_track,
                                        pattern.reordering_power, pattern.reordering_FAR)
                    # AST has already been changed by replace(). However, post 
                    # processing fails. That AST might be wrong. So abort 
                    dprintln(1, 2, "Post-processing failed.")
                    throw(PostPatternReplacementFailure(pattern))
                end
            else
                assert(isa(pattern.post_processing, AbstractPatternAction))
                if !pattern.post_processing.caller(pattern.post_processing, ast, call_sites, pattern.fknob_creator,
                                        pattern.fknob_deletor, pattern.matrices_to_track,
                                        pattern.reordering_power, pattern.reordering_FAR)
                    # AST has already been changed by replace(). However, post 
                    # processing fails. That AST might be wrong. So abort 
                    dprintln(1, 2, "Post-processing failed.")
                    throw(PostPatternReplacementFailure(pattern))
                end
            end
        end

        return true
    end
    return false
end

@doc """
Match and replace AST with the given patterns in top-down fashion, i.e. a parent
tree node in the AST is processed before any child node of it.
"""
function top_down_match_replace_an_expr_pattern(
    ast        :: Expr,
    call_sites :: CallSites,
    patterns   :: Vector #{ExprPattern}
)
    # Match against the Expr. Replacement may happen on the expression.
    for pattern in patterns
        success = match_replace(pattern, ast, call_sites)
        if success
            break
        end
    end
               
    # Now match against each arg. Replacement may happen on the args.
    for arg in ast.args
        if typeof(arg) <: Expr
            top_down_match_replace_an_expr_pattern(arg, call_sites, patterns)
        end
    end            
end

@doc """
Match and replace AST with the given patterns in bottom_up fashion, i.e. a parent
tree node in the AST is processed after all children nodes of it.
"""
function bottom_up_match_replace_an_expr_pattern(
    ast        :: Expr,
    call_sites :: CallSites,
    patterns   :: Vector #{ExprPattern}
)
    # Match against each arg first. Replacement may happen on the args.
    for arg in ast.args
        if typeof(arg) <: Expr
            bottom_up_match_replace_an_expr_pattern(arg, call_sites, patterns)
        end
    end
    
    # Match against the Expr. Replacement may happen on the expression.
    for pattern in patterns
        success = match_replace(pattern, ast, call_sites)
        if success
            break
        end
    end
end

@doc """
Match an expression pattern and do replacement.
"""
function match_replace_an_expr_pattern(
    ast        :: Any,
    call_sites :: CallSites
)    
    if typeof(ast) <: Expr
        # First, do non-splitting patterns. Either top-down or bottom-up
        # matching should be OK. The only difference is that the patterns for
        # each may differ a little. Here we do top-down, for which patterns
        # can be designed based on the shapes of the user source code. In
        # contrast, bottom-up patterns would be based on both the user source
        # code shape and the other patterns that might have been applied.
        if trace_call_replacement
            println("\tTop down matching ast:", ast)
            dsprint(200, 1, call_sites.symbol_info, ast)
        end
        top_down_match_replace_an_expr_pattern(ast, call_sites, 
                                               call_sites.extra.non_splitting_patterns)

        if trace_call_replacement
            println("\tBottom up matching ast:", ast)
            dsprint(200, 1, call_sites.symbol_info, ast)
        end

        # Now do splitting patterns. They should be applied bottom up, because
        # every internal tree node is ensured to be processed. In constrast, 
        # if top-down, a subtree might be hoisted up, and its internal nodes
        # won't be able to be processed subsequently.
        # In addition, we can do non-splitting patterns as well, because some
        # of them are not based source code shapes, but based on other patterns'
        # results. So in effect, we do all patterns once bottom up.
        bottom_up_match_replace_an_expr_pattern(ast, call_sites, 
                                                call_sites.patterns)
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
        if match_skeletons(first_expr,  pattern.first_skeleton,  call_sites, call_sites.extra.live_in_before_prev_expr, first_expr, second_expr) &&
           match_skeletons(second_expr, pattern.second_skeleton, call_sites, call_sites.extra.live_in_before_prev_expr, first_expr, second_expr)
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
            ast_changed        = replace(pattern.new_first_skeleton, first_expr, call_sites, first_expr_backup, second_expr_backup, true)
            ast_changed       |= replace(pattern.new_second_skeleton, second_expr, call_sites, first_expr_backup, second_expr_backup, false)

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
Check if replacing the (symbolic) argument may split the current statement.
So far, only replacing such an argument like :t1!2 can split a statement.
"""
function may_split_statement(
    arg :: Any
)
    if typeof(arg) == Symbol
        arg_string = string(arg)
        if length(arg_string) > 1 && (arg_string[1] == 't')
            return true
        end
    elseif typeof(arg) == Expr
        for a in arg.args
            if may_split_statement(a)
                return true
            end
        end
    end
    return false
end    

@doc """
Separate non-splittable from splittable ExprPatterns. Record the results
in the call sites's extra.
"""
function separate_expr_patterns(
    call_sites :: CallSites
)
    for pattern in call_sites.patterns
        splittable = false
        for arg in pattern.substitute
            if may_split_statement(arg)
                splittable = true
                break
            end
        end
        if !splittable
            push!(call_sites.extra.non_splitting_patterns, pattern)
        end
    end
end

@doc """ 
Pattern match and replace the code in the specified basic blocks that is 
functionally equivalent to SpMV, SpMV!, dot, WAXPBY, WAXPBY!, etc. with calls
to the corresponding SPMP library functions.
"""
function replace_calls(
    bb_indices  :: Set{BasicBlockIndex},
    call_sites  :: CallSites,
    liveness    :: Liveness, 
    cfg         :: CFG
)
    blocks = cfg.basic_blocks
    for bb_idx in bb_indices
        bb                  = blocks[bb_idx]
        call_sites.extra.bb = bb
        prev_expr           = nothing
        statements          = bb.statements
        for stmt_idx in 1 : length(statements)
            call_sites.extra.stmt_idx = stmt_idx
            stmt                      = statements[stmt_idx]
            expr                      = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            if stmt.index < 1
                # This is a new statement generated by, e.g. context-sensitive
                # optimization. It does not have liveness info.  
                call_sites.extra.live_in_before_expr = Set{Sym}()
            else
                call_sites.extra.live_in_before_expr = LivenessAnalysis.live_in(stmt, liveness)
            end
                        
            # Try to pattern match and replace this expression with ExprPatterns.
            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
            
            if prev_expr != nothing
                # Try to merge this and the previous expression
                match_replace_an_two_statements_pattern(prev_expr, expr, call_sites, two_statements_patterns)
            end
            
            prev_expr = expr
            call_sites.extra.live_in_before_prev_expr = call_sites.extra.live_in_before_expr
        end
    end
end

@doc """ 
Pattern match and replace the code that is functionally equivalent to SpMV, SpMV!,
dot, WAXPBY, WAXPBY!, etc. with calls to the corresponding SPMP library functions.
"""
function replace_calls(
    func_region :: FunctionRegion,
    regions     :: Vector{LoopRegion},
    lambda      :: LambdaInfo,
    symbol_info :: Sym2TypeMap,
    liveness    :: Liveness, 
    cfg         :: CFG
)
    dprintln(1, 0, "\nReplacing sparse matrix function calls:")

    if isempty(func_region.symbol_property)
        func_region.symbol_property = find_properties_of_matrices(func_region, lambda, symbol_info, liveness, cfg)
    end
    matrix_properties           = func_region.symbol_property
    call_sites                  = CallSites(Set{CallSite}(), func_region, 
                                            lambda, symbol_info, 
                                            liveness, expr_patterns, 
                                            Vector{Action}(), 
                                            CallReplacementExtra(matrix_properties))

    # Separate non-splittable from splittable patterns.
    separate_expr_patterns(call_sites)
    if !use_splitting_patterns
        call_sites.patterns = call_sites.extra.non_splitting_patterns
    end


    # Separate basic blocks that do not belong to any loop region
    non_loop_bb_indices = Set{BasicBlockIndex}()
    for (bb_idx, bb) in cfg.basic_blocks
        inside_loop = false
        for region in regions
            L = region.loop
            if in(bb_idx, L.members)
                inside_loop = true
                break
            end
        end
        if !inside_loop
            push!(non_loop_bb_indices, bb_idx)
        end
    end
    
    # Replace calls in non loop basic blocks.
    replace_calls(non_loop_bb_indices, call_sites, liveness, cfg)
    
    # Replace calls in each loop region
    for region in regions
        call_sites.region                  = region
        if isempty(region.symbol_property)
            find_properties_of_matrices(region, lambda, symbol_info, liveness, cfg)
        end    
        call_sites.extra.matrix_properties = region.symbol_property 
        replace_calls(region.loop.members, call_sites, liveness, cfg)
    end
    
    return call_sites.actions
end
