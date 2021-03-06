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

@doc """
A function knob is a symbol in compile time. For this symbol, we would like to
record a type for it in the Sym2TypeMap. This abstract type is for this purpose
only.
"""
abstract FunctionKnob

@doc """
The CallSites' extra field for context info discovery.
"""
type ContextExtra
    # Environment before the AST walking
    matrix_properties        :: Symexpr2PropertiesMap
    bb2depth                 :: Dict{BasicBlock, Int}

    # Flag if ast may be under change. We visite the AST twice. It is true
    # the first time (allowing pattern replacement; knobs are generated as well),
    # and false the second time (allowing gathering the knobs of the AST nodes).
    ast_may_change           :: Bool

    # Infomration gathered when AST is walked the first time
    reordering_decider_fknob :: Symbol # Remember fknob instead of AST here as AST may change due to call replacement, etc.
    reordering_decider_power :: Int    # 0 is no power: not to make reordering decision
    reordering_FAR           :: Vector # First Arrays Reordered by the decider
    expr2fknob               :: Dict{Expr, Symbol}
    fknob2expr               :: Dict{Symbol, Expr}
    fknob2mknobs             :: Dict{Symbol, Vector} # Vector[matrix knobs]
    fknob2ctor_dtor          :: Dict{Symbol, Tuple}  # Tuple(fknob_creator, fknob_deletor)
    matrix2mknob             :: Dict{Symexpr, Symbol}
    mknob2matrix             :: Dict{Symbol, Symexpr}
    derivatives              :: Set{Tuple} # Tuple(matrixA, derivative relation, matrixB)

    # The knobs gathered when AST is walked the second time.
    matrix_knobs             :: Set{Symbol}
    function_knobs           :: Set{Symbol}
 
    # Scratch variables.
    bb                       :: Any              # BasicBlock
    stmt_idx                 :: StatementIndex
    prev_stmt_idx            :: StatementIndex
    live_in_before_prev_expr :: Set{Sym}
    live_in_before_expr      :: Set{Sym}

    # Some patterns (those contain :t1!2, etc.) may hoist some subtrees of the
    # current statement before it. That splits one statement into more than one.
    # Such patterns should be matched at the last, because otherwise, other 
    # patterns may not be able to match what they should: they cannot find the
    # subtrees to match, which are no longer in the same statement.
    non_splitting_patterns  :: Vector{Pattern}

    ContextExtra(_matrix_properties, _bb2depth, _ast_may_change) = new(
            _matrix_properties, _bb2depth, _ast_may_change, :nothing, 0, [],
            Dict{Expr, Symbol}(), Dict{Symbol, Expr}(), Dict{Symbol, Vector}(),
            Dict{Symbol, Tuple}(), Dict{Symexpr, Symbol}(),
            Dict{Symbol, Symexpr}(), Set{Tuple}(),
            Set{Symbol}(), Set{Symbol}(), nothing, 0, 0, Set{Sym}(), Set{Sym}(),
            Vector{Pattern}())
end

# Below are the patterns that we would like to replace a function with another,
# which is semantically equivalent, but is context-sensitive
# For short, CS represents "Context Sensitive".

@doc """
Pre-processing function of CS_fwd/bwdTriSolve!_pattern. It checks if the L/U
matrix has a proxy structure available or not. That proxy matrix has to be a
symmetric matrix (not just lower or upper part). This is due to the special 
requirement of the library to fast build a dependence graph. If the matrix is just
lower or upper part, the library building the graph would be twice slower.
"""
function CS_fwdBwdTriSolve!_check(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString
)
    arg2 = ast.args[2]
    if typeof(arg2) == SymbolNode
        L_or_U = arg2.name
    else
        L_or_U = arg2
    end
    assert(typeof(L_or_U) <: Symexpr)
    
    if !haskey(call_sites.extra.matrix_properties, L_or_U) ||
        (call_sites.extra.matrix_properties[L_or_U].lower_of == nothing &&
         call_sites.extra.matrix_properties[L_or_U].upper_of == nothing)
        return false
    end
    
    # It is not possible to be a lower of some matrix, and also an upper of some
    # (other) matrix.
    assert(call_sites.extra.matrix_properties[L_or_U].lower_of == nothing ||
           call_sites.extra.matrix_properties[L_or_U].upper_of == nothing)
    
    # Check it is lower or upper part of a symmetric matrix.
    if call_sites.extra.matrix_properties[L_or_U].lower_of != nothing
        M = get_symexpr(call_sites.extra.matrix_properties[L_or_U].lower_of)
    else
        M = get_symexpr(call_sites.extra.matrix_properties[L_or_U].upper_of)
    end

    if !haskey(call_sites.extra.matrix_properties, M) ||
       !call_sites.extra.matrix_properties[M].constant_structured
       !call_sites.extra.matrix_properties[M].is_structure_symmetric
        return false
    end


    if typeof(L_or_U) == Symbol && typeof(M) == Symbol
        push!(call_sites.extra.derivatives, (L_or_U, DERIVATIVE_TYPE_SYMMETRIC, M))
    end

    return true
end

@doc """ 
Post-processing function of a pattern: Gather context sensitive info of the call
site at the given AST.
"""
function gather_context_sensitive_info(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: AbstractString,
    fknob_deletor     :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    assert(call_sites.extra.ast_may_change)

    # Create a function-specific knob. 
    fknob                                   = new_symbol("fknob")
    call_sites.extra.expr2fknob[ast]        = fknob
    call_sites.extra.fknob2expr[fknob]      = ast
    call_sites.extra.fknob2ctor_dtor[fknob] = (fknob_creator, fknob_deletor)
    call_sites.extra.fknob2mknobs[fknob]    = []
    
    # Record the new symbol with a FunctionKnob type so that we may distinguish
    # it from a normal argument. 
    symbol_info        = call_sites.symbol_info
    symbol_info[fknob] = FunctionKnob

    # Get the previous and the current expression that would be used for
    # argument replacement below.
    bb            = call_sites.extra.bb
    stmt_idx      = call_sites.extra.stmt_idx
    prev_stmt_idx = call_sites.extra.prev_stmt_idx
    cur_expr      = bb.statements[stmt_idx].expr
    prev_expr     = prev_stmt_idx == 0 ? nothing : bb.statements[prev_stmt_idx].expr

    # Create matrix-specific knobs for the matrices inputs.
    args = ast.args
    for arg in matrices_to_track
        if arg == :r
            M = ast
        else
            new_arg = replacement_arg(arg, args, ast, call_sites, prev_expr, cur_expr)
            assert(typeof(new_arg) == Symbol || typeof(new_arg) == SymbolNode ||
                   typeof(new_arg) == GenSym)
            M = ((typeof(new_arg) == SymbolNode) ? new_arg.name : new_arg)
        end
        if !haskey(call_sites.extra.matrix2mknob, M)
            mknob = new_symbol(string("mknob", typeof(M) == Expr ? "Expr" : string(M)))
            call_sites.extra.matrix2mknob[M]     = mknob
            call_sites.extra.mknob2matrix[mknob] = M
        end
        push!(call_sites.extra.fknob2mknobs[fknob], call_sites.extra.matrix2mknob[M])
    end

    bb     = call_sites.extra.bb
    depth  = call_sites.extra.bb2depth[bb]
    weight = 2 << depth 
    if reorder_enabled && depth > 0 && call_sites.extra.reordering_decider_power < (reordering_power * weight)
        # The decider should be in a basic block that is an immediate member of
        # a loop region.
        # The decider should not be in the function region: we want one
        # reordering of data to impact many loop iterations, not sequential code.
        call_sites.extra.reordering_decider_fknob = fknob
        call_sites.extra.reordering_decider_power = reordering_power * weight
        call_sites.extra.reordering_FAR           = []
        for arg in reordering_FAR
            new_arg = replacement_arg(arg, args, ast, call_sites, prev_expr, cur_expr)
            assert(typeof(new_arg) == SymbolNode || typeof(new_arg) <: Symexpr)
            M = ((typeof(new_arg) == SymbolNode) ? new_arg.name : new_arg)
            push!(call_sites.extra.reordering_FAR, M)
        end
    end

    return true
end

@doc """ 
Pre-processing function of CS_ADAT_pattern: check it is A*D*A'; D is a diagonal
array; and A has constant value. 
"""
function CS_ADAT_check(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString
)
    A = ast.args[2]
    assert(typeof(A) == SymbolNode)
    assert(typeof(ast.args[4].args[2]) == SymbolNode)
    if A.name != ast.args[4].args[2].name
        return false
    end

    if !haskey(call_sites.extra.matrix_properties, A.name) ||
        !call_sites.extra.matrix_properties[A.name].constant_valued
        return false
    end

    # TODO: check D is diagonal after structure analysis is done.
    return true
end

@doc """ 
Post-processing function of CS_ADAT_pattern.
This will create an action to transpose A before the loop head. Usually we do
not create any action during this phase (the first walk of the AST). This is
an exception.
"""
function CS_ADAT_post_replacement(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: AbstractString,
    fknob_deletor     :: AbstractString,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    # We need to replace A * D * A' into ADB(A, D, A', fknob). At this 
    # moment, it is in the middle form of ADB(A, D, A). That means:
    # (1) Make a symbol AT. Insert AT = A' before the loop (This is to hoist A'
    #     out of loop)
    # (2) Replace ADB(A, D, A) as ADB(A, D, AT).
    # fknob is not added for now, which will be added automatically later.
    action = (typeof(call_sites.region) == FunctionRegion) ? 
        InsertBeforeOrAfterStatement(Vector{Statement}(), true, call_sites.region.entry, 1) :
        InsertBeforeLoopHead(Vector{Statement}(), call_sites.region.loop, true)
    push!(call_sites.actions, action)

    A = ast.args[4]
    assert(typeof(A) == SymbolNode)
    properties = call_sites.extra.matrix_properties[A.name]
    AT = Symbol(string("__", string(A.name), "T__"))
    stmt = Statement(0, Expr(:(=), AT, Expr(:call, GlobalRef(Main, :ctranspose), A)))
    push!(action.new_stmts, stmt)

    ast.args[4] = SymbolNode(AT, A.typ)

    # AT has the same properties as A.
    call_sites.extra.matrix_properties[AT]               = properties
    call_sites.extra.matrix_properties[AT].is_single_def = false # AT definition is out of the region, and not counted as def.

    push!(call_sites.extra.derivatives, (A.name, DERIVATIVE_TYPE_TRANSPOSE, AT))

    # The result is a symmetric matrix.
    # ISSUE: in analysis like constant value/structure analysis, should we 
    # consider each AST node instead of only symbols? If that is the case, we
    # do not have to make properties here.
    call_sites.extra.matrix_properties[ast]                        = MatrixProperties()
    call_sites.extra.matrix_properties[ast].constant_valued        = false #TODO: should determined by if A and D are constant valued
    call_sites.extra.matrix_properties[ast].constant_structured    = true
    call_sites.extra.matrix_properties[ast].is_symmetric           = false #true
    call_sites.extra.matrix_properties[ast].is_structure_symmetric = false #true
    call_sites.extra.matrix_properties[ast].is_structure_only      = false
    call_sites.extra.matrix_properties[ast].is_single_def          = true 

    return gather_context_sensitive_info(ast, call_sites, fknob_creator, 
                                         fknob_deletor, matrices_to_track,
                                         reordering_power, reordering_FAR)
end

@doc """ 
Pre-processing function of SpSquareWithEps_pattern1: check it is X*X',
and X is symmetric valued. 
"""
function SpSquareWithEps_check(
    ast           :: Expr,
    call_sites    :: CallSites,
    fknob_creator :: AbstractString,
    fknob_deletor :: AbstractString
)
    X = ast.args[2]
    assert(typeof(X) == SymbolNode)
    assert(typeof(ast.args[3].args[2]) == SymbolNode)
    if X.name != ast.args[3].args[2].name
        return false
    end

    if !haskey(call_sites.extra.matrix_properties, X.name) ||
        !call_sites.extra.matrix_properties[X.name].is_symmetric
        return false
    end

    return true
end

@doc """
The following patterns are to match the following source code
    for 
        ...
        B = A*D*A' // A and D::SparseMatrixCSC{Float64,Int32}, D is diagnoal matrix
        ...
        R = cholfact_int32(B)
        ...
        dy = R\t2

Then for each function call site, create a fknob. For each of them, for each
SparseMatrix input and output, add a mknob to represent them. Even though
a matrix might have not been allocated yet, we can still represent it
with its mknob. The purpose of the fknobs are to let us find the mkobs shared
accross functions. 

Then for each function call site, we check some properties of the input mknobs;
if they hold, we check th existence of some data of the output mknob; if they do
not exist, create them; otherwise, directly use them.

This procedure is general, and not limited to this specific case.

The above code is transformed into 
    # dy is used in cholmod_factor_inverse_divide(), and the pattern knows that
    # we need to create dy ahead of time, if dy is NOT live befor the loop.
    dy = Array(Cdouble, m)

    create mknob for A, B, R 
    create fknob for ADB, cholfact_int32, and cholmod_factor_inverse_divide
    add mknob_A, B into fknob_ADB
    add mknob_B, R into fknob_cholfact_int32
    add mknob_R into fknob_cholmod_factor_inverse_divide
    for
        ...
        # Sparso.ADB(A', D, A, fknob_ADB)
        #   if fknob_ADB == NULL
        #       return A * D * A'
        #   get the input mknob_A from fknob_ADB
        #   if !mknob_A->constant_structured (This check is actually already done by pattern matching)
        #       return A * D * A'
        #   get the output mknob_B from fknob_ADB
        #   if mknob_B->structure_proxy == NULL
        #       ADAT = adb_inspect(A', A), 
        #       mknob_B->structure_proxy = ADAT
        #       mknob_B->matrix = ADAT
        #   d = diag(D)
        #   adb_execute!(mknob_B->matrix, A', A, d)
        B = Sparso.ADB(A', D, A, fknob_ADB)
        ...
        # Sparso.cholfact_int32(B, fknobcholfact_int32)
        #   if fknobcholfact_int32 == NULL
        #       return cholfact_int32(B)
        #   get the input mknob_B from fknob_cholfact_int32
        #   if mknob_B->structure_proxy == NULL || mknob_B->matrix == NULL
        #       return cholfact_int32(B)
        #   get the output mknob_R from fknob_cholfact_int32
        #   if mknob_R->dss_handle == NULL
        #       dss_handle = dss_analyze(mknob_B->structure_proxy)
        #       mknob_R->dss_handle = dss_handle
        #   dss_factor(mknob_R->dss_handle, mknob_B->matrix)
        #   BUT HOW TO RETURN R?
        
        R = Sparso.cholfact_int32(B, fknob_cholfact_int32)
        ...
        
        # cholmod_factor_inverse_divide(dy, R, t2, fknob_cholmod_factor_inverse_divide)
        #   if fknob_cholmod_factor_inverse_divide == NULL
        #       return R \ t2
        #   get the input mknob_R from fknob_cholmod_factor_inverse_divide
        #   if mknob_R->dss_handle == NULL
        #       return R \ t2
        #   else
        #       opt = MKL_DSS_DEFAULTS
        #       dss_solve!(mknob_R->dss_handle, t2, dy)

        dy = cholmod_factor_inverse_divide(R, t2, fknob_cholmod_factor_inverse_divide)
"""

const CS_transpose_pattern = ExprPattern(
    "CS_transpose_pattern",
    (:call, GlobalRef(Main, :ctranspose), SparseMatrixCSC),
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

const CS_ADAT_pattern = ExprPattern(
    "CS_ADAT_pattern",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC),
    (nothing, nothing, nothing, nothing, CS_transpose_pattern),
    CS_ADAT_check,
    (:call, GlobalRef(Sparso, :ADB), :a2, :a3, :a2),
    CS_ADAT_post_replacement,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:r, :a2, :a3, :a4),
    0,
    ()
)

const CS_cholfact_int32_pattern = ExprPattern(
    "CS_cholfact_int32_pattern",
    (:call, GlobalRef(Main, :cholfact_int32), SparseMatrixCSC{Float64, Int32}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :cholfact_int32), :a2),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a2,),
    0,
    ()
)

const CS_cholsolve_pattern = ExprPattern(
    "CS_cholsolve_pattern",
    (:call, GlobalRef(Main, :\), Factorization{Float64}, Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :cholfact_inverse_divide), :a2, :a3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a2,),
    0,
    ()
)


@doc """ 
Pattern for foward triangular solver.
"""
const CS_fwdTriSolve!_pattern = ExprPattern(
    "CS_fwdTriSolve!_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), GlobalRef(Base, :SparseMatrix), QuoteNode(:fwdTriSolve!)),
      SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    CS_fwdBwdTriSolve!_check,
    (:call, GlobalRef(Sparso, :fwdTriSolve!), :a2, :a3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a2,),
    100,
    (:a2, :a3)
)

@doc """ 
Pattern for forward triangular solver (A lower triangular matrix \ a vector). The
matrix's structure is the lower part of a constant symmetric matrix's structure.
"""
const CS_fwdTriSolve!_backslash_pattern = ExprPattern(
    "CS_fwdTriSolve!_backslash_pattern",
    (:(=), Vector, 
      Expr(:call, GlobalRef(Main, :\), 
            AD(SparseMatrixCSC, (SA_LOWER_OF, SA_CONST_STRUCTURED | SA_SYMM_STRUCTURED)),
            Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :fwdTriSolve!), :a1, :a2_2, :a2_3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    100,
    (:a3, :a2)
)
 
@doc """
Pattern for backward triangular solver.
"""
const CS_bwdTriSolve!_pattern = ExprPattern(
    "CS_bwdTriSolve!_pattern",
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), GlobalRef(Base, :SparseMatrix), QuoteNode(:bwdTriSolve!)),
      Any, Vector), #SparseMatrixCSC, Vector), # Somehow, type inference does not figure out U's type is SparseMatrixCSC
    (:NO_SUB_PATTERNS,),
    CS_fwdBwdTriSolve!_check,
    (:call, GlobalRef(Sparso, :bwdTriSolve!), :a2, :a3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a2,),
    90,
    (:a2, :a3)
)

@doc """ 
Pattern for backward triangular solver (A upper triangular matrix \ a vector). The
matrix's structure is the upper part of a constant symmetric matrix's structure.
"""
const CS_bwdTriSolve!_backslash_pattern = ExprPattern(
    "CS_bwdTriSolve!_backslash_pattern",
    (:(=), Vector, 
      Expr(:call, GlobalRef(Main, :\), 
            AD(SparseMatrixCSC, (SA_UPPER_OF, SA_CONST_STRUCTURED | SA_SYMM_STRUCTURED)),
            Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :bwdTriSolve!), :a1, :a2_2, :a2_3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    100,
    (:a3, :a2)
)

@doc """
The following patterns are for SpMV. The difference from the SpMV patterns in
call-replacement.jl is that here we will add mknobs and fknob to SpMV.
"""

@doc """ A * x => SpMV(1, A, x) """
const CS_SpMV_pattern1 = ExprPattern(
    "CS_SpMV_pattern1",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV), 1, :a2, :a3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    40,
    (:a3, :a4, :r)
)

@doc """ a * A * x + g => SpMV(a, A, x, 0, x, g) """
const CS_SpMV_pattern2 = ExprPattern(
    "CS_SpMV_pattern2",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, number_times_matrix_vector_pattern, nothing),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV), :a2_2, :a2_3, :a2_4, 0, :a2_4, :a3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    40,
    (:a3, :a4, :r)
)

@doc """ a * A * x => SpMV(a, A, x) """
const CS_SpMV_pattern3 = ExprPattern(
    "CS_SpMV_pattern3",
    (:call, GlobalRef(Main, :*), Number, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV), :a2, :a3, :a4),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    40,
    (:a3, :a4, :r)
)

@doc """
SpMV(a, A, x) + g => SpMV(a, A, x, 0, x, g)
Note: its reordering power is bigger than CS_SpMV_pattern3, because
CS_SpMV_pattern3 generates a SpMV_3_parameters_pattern, which is a sub-pattern
of this pattern. Making reordering power of this pattern bigger ensures that 
when the sub-pattern is replaced, the reordering decider becomes the AST
generated by this pattern, and no more being that by the sub-pattern.
"""
const CS_SpMV_pattern4 = ExprPattern(
    "CS_SpMV_pattern4",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, SpMV_3_parameters_pattern, nothing),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV),
     :a2_2, :a2_3, :a2_4, 0, :a2_4, :a3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    50,
    (:a3, :a4, :r)
)

@doc """
SpMV(a, A, x) + b => SpMV(a, A, x, 1, b, 0)
Note: its reordering power is bigger than CS_SpMV_pattern3, because
CS_SpMV_pattern3 generates a SpMV_3_parameters_pattern, which is a sub-pattern
of this pattern. Making reordering power of this pattern bigger ensures that 
when the sub-pattern is replaced, the reordering decider becomes the AST
generated by this pattern, and no more being that by the sub-pattern.
"""
const CS_SpMV_pattern5 = ExprPattern(
    "CS_SpMV_pattern5",
    (:call, GlobalRef(Main, :+), Vector, Vector),
    (nothing, nothing, SpMV_3_parameters_pattern, nothing),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV),
     :a2_2, :a2_3, :a2_4, 1, :a3, 0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    50,
    (:a3, :a4, :a6, :r)
)

@doc """
SpMV(a, A, x) - b => SpMV(a, A, x, -1, b, 0)
Note: its reordering power is bigger than CS_SpMV_pattern3, because
CS_SpMV_pattern3 generates a SpMV_3_parameters_pattern, which is a sub-pattern
of this pattern. Making reordering power of this pattern bigger ensures that 
when the sub-pattern is replaced, the reordering decider becomes the AST
generated by this pattern, and no more being that by the sub-pattern.
"""
const CS_SpMV_pattern6 = ExprPattern(
    "CS_SpMV_pattern6",
    (:call, GlobalRef(Main, :-), Vector, Vector),
    (nothing, nothing, SpMV_3_parameters_pattern, nothing),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV),
     :a2_2, :a2_3, :a2_4, -1, :a3, 0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    50,
    (:a3, :a4, :a6, :r)
)

@doc """
    (-1/a)*A*x + b*y => SpMV(-1/a, A, x, b, y)
"""
const CS_SpMV_pattern7 = ExprPattern(
    "CS_SpMV_pattern7",
    (:call, GlobalRef(Main, :+), 
      Expr(:call, GlobalRef(Main, :/), 
        Expr(:call, GlobalRef(Main, :-),      
            Expr(:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector)),
        Number),     
      Expr(:call, GlobalRef(Main, :*), Number, Vector)), 
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV),
      TypedExprNode(Float64, :call, GlobalRef(Main, :/), -1, :a2_3), :a2_2_2_2, :a2_2_2_3, :a3_2, :a3_3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    70,
    (:a3, :a4, :a6, :r)
)            

@doc """
The following patterns are potentially related with SpMV!. In SpMP library, SpMV!
is implemented in such a way that the result and the inputs are written and read
in parallel. Thus the result and the inputs cannot share space to avoid any race
condition.

There are 3 cases:

(Case 1): y = A *x ==> y = SpMV(A, x)

This is equivalent to
        result = allocate memory of size(A, 1)
        SpMV!(result, A, x)
        y = result
The first two statements have been realized by SpMV. Thus we can simply generate
code like this
        y = SpMV(A, x)

One might think it is great to hoist the allocation out of the current region (
a loop). However, that can cause subtle issues. Consider a simplified pagerank:
        y = A * x
        if |y - x| < err exit
        x = y

If we hoist the allocation out, the code becomes:
    Before the loop:
        result = allocate memory of size(A, 1)
    In the loop
        SpMV!(result, A, x)
        y = result
        if |y - x| < err exit
        x = y
After the first iteration, x, y, and result share the same memory. Thus |y-x|
is always 0 afterwards.

(Case 2) x = A * x  ==> x = SpMV(A, x) or SpMV!(x, A, x) with library extension.

Here we are sure that x has been allocated memory. There are three ways we can 
generate code: 

CodeGen 21: generate 
        x = SpMV(A, x)
The overhead is: memory allocation and copy every time SpMV! is called.

CodeGen 22: generate 
        SpMV!(x, A, x)
This requires that the library is extended to allow that the result and 
the input are aliased. To do that, library identifies such an alias, and 
creates a temporary memory space t as result of A * x, and run SpMV!(t, A, x),
where it can write and read inputs in parallel, and finally copy t's content to
x. The overhead is the same: memory allocation and copy every time SpMV! is
called.

CodeGen 23: generate
    Before the loop:
        result = allocate memory of size(A, 1)
        assert(result and x have different spaces)
    In the loop:
        SpMV!(result, A, x)
        tmp = x
        x = result
        result = tmp
        assert(result and x have different spaces)
This has 1 allocation, and no copy back.

So far, we do not have Case 2 in any benchmark yet. So let us leave it for future.
Simply follow CodeGen 21 or 22 for now.

(Case 3) A_mul_B!(y, A, x) ==> SpMV!(y, A, x) with library extension.

This is equivalent to y = A * x, but y is sure to have been allocated memory
before.

If y is lexically x, it is a special case A_mul_B!(x, A, x). Same as Case 2.

Otherwise, we need to know if y and x are aliased (share the same memory). 
If they are aliased, that is Case 2.
If they are not, simply call SpMV!(y, A, x) as usual.

So far, we do not rely on static alias analysis. So probably we should extend
libary to identify and allow the result and the input being aliased. See notes 
for CodeGen 22.
"""

@doc """
A_mul_B!(y, A, x) = SpMV!(y, A, x)
"""
const CS_SpMV!_pattern1 = ExprPattern(
    "CS_SpMV!_pattern1",
    (:call, GlobalRef(Main, :A_mul_B!), Vector, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :a2, :a3, :a4),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a3,),
    50,
    (:a3, :a4, :a2) # Put seed (A) as the first
)

@doc """ z = SpMV(a, A, x, b, y, g), z is x or y => SpMV!(z, a, A, x, b, y, g) """
const CS_SpMV!_pattern2 = ExprPattern(
    "CS_SpMV!_pattern2",
    (:(=), Vector, Vector),
    (nothing, nothing, SpMV_6_parameters_pattern),
    # NOTE: this is a HACK! We should not do this without knowing that z has memory
    # allocated. This is only to enable using SpMV! for performance.
    # TODO: uncomment this condition below after Linxiang's analysis works
    # LHS_in_RHS,
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :a1, :a2_2, :a2_3, :a2_4, :a2_5, :a2_6, :a2_7),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a4,),
    80,
    (:a4, :a2, :a5, :a7) # Put seed (A) at the first
)

@doc """ w = a * A * x => SpMV!(w, a, A, x, 0, x, 0) """
const CS_SpMV!_pattern3 = ExprPattern(
    "CS_SpMV!_pattern3",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY),
        Expr(:call, GlobalRef(Sparso, :SpMV),
              Number, SparseMatrixCSC, Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :a1, :a2_2, :a2_3, :a2_4, 0.0, :a1, 0.0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a4,),
    60,
    (:a4, :a2, :a5) # Put seed (A) at the first
)            

@doc """ A_mul_B!(a, A, x, b, y) => SpMV!(y, a, A, x, b, y , 0) """
const CS_SpMV!_pattern4 = ExprPattern(
    "CS_SpMV!_pattern4",
    (:call, GlobalRef(Main, :A_mul_B!), Number, SparseMatrixCSC, Vector, Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :a6, :a2, :a3, :a4, :a5, :a6, 0.0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a4,),
    70,
    (:a4, :a2, :a5, :a7) # Put seed (A) at the first
)

@doc """ z = SpMV(a, A, x, b, y), z has dedicated memory => SpMV!(z, a, A, x, b, y) """
const CS_SpMV!_pattern5 = ExprPattern(
    "CS_SpMV!_pattern5",
    (:(=), AD(Vector, SA_HAS_FREE_MEMORY),
        Expr(:call, GlobalRef(Sparso, :SpMV),
             Number, SparseMatrixCSC, Vector, Number, Vector)),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :a1, :a2_2, :a2_3, :a2_4, :a2_5, :a2_6, 0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a4,),
    80,
    (:a4, :a2, :a5, :a7) # Put seed (A) at the first
)

@doc """ spmatmul_witheps(X, X', eps), X is symmetric => SpSquareWithEps(X, eps) """
const SpSquareWithEps_pattern1 = ExprPattern(
    "SpSquareWithEps_pattern1",
    (:call,   GlobalRef(Main, :spmatmul_witheps), SparseMatrixCSC, SparseMatrixCSC,      Number),
    (nothing, nothing,                            nothing,         CS_transpose_pattern, nothing),
    SpSquareWithEps_check,
    (:call, GlobalRef(Sparso, :SpSquareWithEps),
     :a2, :a4),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:a2,),
    0,
    ()
)

@doc """" Patterns that will actually transform the code. """
CS_transformation_patterns = [
    CS_ADAT_pattern,
    CS_cholfact_int32_pattern,
    CS_cholsolve_pattern,
    CS_fwdTriSolve!_pattern,
    CS_fwdTriSolve!_backslash_pattern,
    CS_bwdTriSolve!_pattern,
    CS_bwdTriSolve!_backslash_pattern,
    CS_SpMV_pattern7,
    CS_SpMV_pattern1,
    CS_SpMV_pattern2,
    CS_SpMV_pattern3,
    CS_SpMV_pattern4,
    CS_SpMV_pattern5,
    CS_SpMV_pattern6,
    CS_SpMV!_pattern1,
    CS_SpMV!_pattern2,
    CS_SpMV!_pattern3,
    # CS_SpMV!_pattern4] # Do not handle y = A * x + y case for now
    CS_SpMV!_pattern5,
    SpSquareWithEps_pattern1
]

@doc """
    Sparso.SpMV!(z, a, A, x, b, y, r)
    z = z .* u
=>
    Sparso.SpMV!(z, a, A, x, b, y, r, u)
    z = z # Useless
"""
const CS_SpMV!_two_statements_pattern1 = TwoStatementsPattern(
    "CS_SpMV!_two_statements_pattern1",
    (:call, GlobalRef(Sparso, :SpMV!),
      Vector, Number, SparseMatrixCSC, Vector, Number, Vector, Number),
    (:(=), :f2, Expr(:call, GlobalRef(Main, :.*), :f2, Vector)),
    do_nothing,
    (:call, GlobalRef(Sparso, :SpMV!),
     :f2, :f3, :f4, :f5, :f6, :f7, :f8, :s2_3),
    (:(=), :f2, :f2),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:f4,),
    100,
    (:f4, :f2, :f5, :f7, :f9) # Put seed (A) at the first
)

CS_two_statements_patterns = [
    CS_SpMV!_two_statements_pattern1,
]


@doc """
Create statements that will create a matrix knob for matrix M.
"""
function create_new_matrix_knob(
    new_stmts  :: Vector{Statement},
    mknob      :: Symbol,
    M          :: Symexpr,
    call_sites :: CallSites
)
    properties = call_sites.extra.matrix_properties[M]
    new_stmt = Expr(:(=), mknob,
                    Expr(:call, GlobalRef(Sparso, :new_matrix_knob), 
                         typeof(M) == Expr ? QuoteNode(symbol("Expr")) : QuoteNode(M),
                         properties.constant_valued, 
                         properties.constant_structured, 
                         properties.is_symmetric,
                         properties.is_structure_symmetric,
                         properties.is_structure_only,
                         properties.is_single_def))
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create a statement to propagate matrix info from one matrix knob to another.
"""
function create_propagate_matrix_info(
    to_mknob   :: Symbol,
    from_mknob :: Symbol,
    bb         :: BasicBlock,
    stmt_idx   :: StatementIndex,
    call_sites :: CallSites
)
    # Insert the propagation after the (assignment) statement.
    action = InsertBeforeOrAfterStatement(Vector{Statement}(), false, bb, stmt_idx)
    stmt = Statement(0, 
                     Expr(:call, GlobalRef(Sparso, :propagate_matrix_info), 
                           to_mknob, from_mknob))
    push!(action.new_stmts, stmt)
    push!(call_sites.actions, action)
end

@doc """
Create statements that will create a function knob for the call site, and add
the function knob to the call as a parameter.
"""
function create_new_function_knob(
    new_stmts     :: Vector{Statement},
    fknob         :: Symbol,
    fknob_creator :: AbstractString,
)
    assert(fknob_creator != "")

    if fknob_creator == "NewFunctionKnob"
        # Create a function knob. It has no private info specific to that function.
        # Call the parameterless version of new_function_knob for cleaner code.
        new_stmt = Expr(:(=), fknob, 
                        Expr(:call, GlobalRef(Sparso, :new_function_knob))
                   )
    else
        # So far, we do not create any special knob. So this case is never reached.
        assert(false)
        new_stmt = Expr(:(=), fknob, 
                        Expr(:call, GlobalRef(Sparso, :new_function_knob),
                              fknob_creator
                        )
                   )
    end
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create statements that add the matrix knob to the function knob.
"""
function create_add_mknob_to_fknob(
    new_stmts :: Vector{Statement},
    mknob     :: Symbol,
    fknob     :: Symbol
)
    new_stmt = Expr(:call, GlobalRef(Sparso, :add_mknob_to_fknob), mknob, fknob)
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create statements that will delete the function knob.
"""
function create_delete_function_knob(
    new_stmts     :: Vector{Statement},
    fknob         :: Symbol,
    fknob_deletor :: AbstractString
)
    assert(fknob_deletor != "")

    if fknob_deletor == "DeleteFunctionKnob"
        # Delete a function knob, which has no private info specific to that function.
        # Call the 1-parameter version of delete_function_knob for cleaner code.
        new_stmt = Expr(:call, GlobalRef(Sparso, :delete_function_knob), 
                         fknob)
    else
        # So far, we have not created any special knob. So this case is never reached.
        assert(false)
        new_stmt = Expr(:call, GlobalRef(Sparso, :delete_function_knob),
                         fknob_deletor, fknob)
    end
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create statements that will delete the matrix knob.
"""
function create_delete_matrix_knob(
    new_stmts :: Vector{Statement},
    mknob     :: Symbol
)
    new_stmt = Expr(:call, GlobalRef(Sparso, :delete_matrix_knob), mknob)
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create statements that set mknob2 as a derivative of mknob1.
"""
function create_set_derivative(
    new_stmts :: Vector{Statement},
    mknob1    :: Symbol,
    relation  :: Int,
    mknob2    :: Symbol
)
    new_stmt = Expr(:call, GlobalRef(Sparso, :set_derivative),
                     mknob1, int2derivative_map[relation], mknob2)
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Visit each expression the loop region. If recursive is true, visit each AST node
of the expression and handle each with the handler. Otherwise, handle the
expression with the handler.
"""
function visit_expressions(
    region                 :: Region,
    cfg                    :: CFG,
    call_sites             :: CallSites,
    recursive              :: Bool,
    handler                :: Function,
    two_statements_handler :: Function = do_nothing
)
    members     = (typeof(region) == FunctionRegion) ? region.members : region.loop.members
    blocks      = cfg.basic_blocks
    symbol_info = call_sites.symbol_info
    liveness    = call_sites.liveness
    for bb_idx in members
        bb                             = blocks[bb_idx]
        statements                     = bb.statements
        call_sites.extra.bb            = bb
        prev_expr                      = nothing
        call_sites.extra.prev_stmt_idx = 0
        for stmt_idx in 1 : length(statements)
            call_sites.extra.stmt_idx = stmt_idx
            stmt                      = statements[stmt_idx]
            expr                      = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            call_sites.extra.live_in_before_expr = LivenessAnalysis.live_in(stmt, liveness)

            if recursive
                CompilerTools.AstWalker.AstWalk(expr, handler, call_sites)
            else
                handler(expr, call_sites)
            end
                
            if prev_expr != nothing && two_statements_handler != do_nothing
                # Try to merge this and the previous expression
                two_statements_handler(prev_expr, expr, call_sites, CS_two_statements_patterns)
            end
            prev_expr                                 = expr
            call_sites.extra.prev_stmt_idx            = stmt_idx
           call_sites.extra.live_in_before_prev_expr = call_sites.extra.live_in_before_expr            
        end
    end
end

@doc """"
Add the function and matrix knobs of the AST node to 
call_sites.extra.matrix_knobs and function_knobs.
"""
function gather_knobs(
    ast        :: Any,
    call_sites :: CallSites
)
    assert(!call_sites.extra.ast_may_change)

    if typeof(ast) != Expr
        return
    end

    if haskey(call_sites.extra.expr2fknob, ast)
        fknob = call_sites.extra.expr2fknob[ast]
        push!(call_sites.extra.function_knobs, fknob)
        
        mknobs = call_sites.extra.fknob2mknobs[fknob]
        union!(call_sites.extra.matrix_knobs, mknobs)

        site = CallSite(ast)
        push!(call_sites.sites, site)
    end
    return nothing
end

@doc """ A handler to visit_expressions(). """
gather_knobs(ast, call_sites :: CallSites, top_level_number, is_top_level, read) =
    gather_knobs(ast, call_sites)

@doc """ Create statements that will delete matrix and function knobs. """
function delete_knobs(
    action         :: Action,
    matrix_knobs   :: Set{Symbol},
    function_knobs :: Set{Symbol},
    call_sites     :: CallSites
)
    for mknob in matrix_knobs
        create_delete_matrix_knob(action.new_stmts, mknob)
    end

    for fknob in function_knobs
        fknob_creator, fknob_deletor = call_sites.extra.fknob2ctor_dtor[fknob]
        create_delete_function_knob(action.new_stmts, fknob, fknob_deletor)
    end
end

@doc """
Create statements that will create and intialize a knob for a matrix/function
before the loop, and statements that will delete all the matrix/function
knobs at each exit of the loop.
"""
function generate_and_delete_knobs(
    region     :: Region,
    call_sites :: CallSites
)
    region         = call_sites.region
    matrix_knobs   = call_sites.extra.matrix_knobs
    function_knobs = call_sites.extra.function_knobs
    derivatives    = call_sites.extra.derivatives
    symbol_info    = call_sites.symbol_info

    action_before_region = (typeof(region) == FunctionRegion) ? 
        InsertBeforeOrAfterStatement(Vector{Statement}(), true, region.entry, 1) :
        InsertBeforeLoopHead(Vector{Statement}(), region.loop, true)
    push!(call_sites.actions, action_before_region)

    # Add statements that generates mknobs and fknobs
    for mknob in matrix_knobs
        M = call_sites.extra.mknob2matrix[mknob]
        create_new_matrix_knob(action_before_region.new_stmts, mknob, M, call_sites)
    end

    for fknob in function_knobs
        fknob_creator, fknob_deletor = call_sites.extra.fknob2ctor_dtor[fknob]
        create_new_function_knob(action_before_region.new_stmts, fknob, fknob_creator)
    end

    # Derivative information is already included in matrix_knob, so we copy it into derivatives
    for mknob in matrix_knobs
        M = call_sites.extra.mknob2matrix[mknob]
        prop = call_sites.extra.matrix_properties[M]
        if prop.lower_of != nothing && haskey(call_sites.extra.matrix2mknob, prop.lower_of)
            assert(prop.upper_of == nothing)
            push!(derivatives, (M, DERIVATIVE_TYPE_LOWER_TRIANGULAR, prop.lower_of))

            # ISSUE: is the following condition right for setting symmetric derivative? Do we
            # need both prop1.is_symmetric && prop1.is_structure_symmetric? Is this derivative
            # symmetric in structure or value?
            # ISSUE: without setting the symmetric derivative, pcg-illustration causes spmp
            # segmentation fault
            prop1 = call_sites.extra.matrix_properties[prop.lower_of]
            if prop1.is_symmetric && prop1.is_structure_symmetric 
                push!(derivatives, (M, DERIVATIVE_TYPE_SYMMETRIC, prop.lower_of))
            end
        elseif prop.upper_of != nothing && haskey(call_sites.extra.matrix2mknob, prop.upper_of)
            push!(derivatives, (M, DERIVATIVE_TYPE_UPPER_TRIANGULAR, prop.upper_of))
            
            # ISSUE: is the following condition right for setting symmetric derivative? Do we
            # need both prop1.is_symmetric && prop1.is_structure_symmetric? Is this derivative
            # symmetric in structure or value?
            # ISSUE: without setting the symmetric derivative, pcg-illustration causes spmp
            # segmentation fault
            prop1 = call_sites.extra.matrix_properties[prop.upper_of]
            if prop1.is_symmetric && prop1.is_structure_symmetric 
                push!(derivatives, (M, DERIVATIVE_TYPE_SYMMETRIC, prop.upper_of))
            end
        end
        if prop.transpose_of != nothing && haskey(call_sites.extra.matrix2mknob, prop.transpose_of)
            push!(derivatives, (M, DERIVATIVE_TYPE_TRANSPOSE, prop.transpose_of))
        end
    end

    # Add statements that define the relationship between matrix knobs
    for (M1, relation, M2) in derivatives
        mknob1 = call_sites.extra.matrix2mknob[M1]
        mknob2 = call_sites.extra.matrix2mknob[M2]
        create_set_derivative(action_before_region.new_stmts, mknob1, relation, mknob2)
    end

    if (typeof(region) == FunctionRegion)
        for pred_bb in region.exit.preds
            len = length(pred_bb.statements)
            assert(len > 0)
            last_stmt = pred_bb.statements[len]
            last_expr = last_stmt.expr
            assert(typeof(last_expr) <: Expr && last_expr.head == :return)

            action  = InsertBeforeOrAfterStatement(Vector{Statement}(), true, pred_bb, len)
            push!(call_sites.actions, action)
            delete_knobs(action, matrix_knobs, function_knobs, call_sites)
        end
    else
        # Create statements that will delete the knobs at region exits
        for exit in region.exits
            action  = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
            push!(call_sites.actions, action)
            delete_knobs(action, matrix_knobs, function_knobs, call_sites)
        end
    end
end

@doc """"
Create statements before the region that add matrix knobs into the function knob.
"""
function add_matrix_knobs_to_function_knobs(
    call_sites :: CallSites
)
    assert(!call_sites.extra.ast_may_change)

    region               = call_sites.region
    action_before_region = (typeof(region) == FunctionRegion) ? 
        InsertBeforeOrAfterStatement(Vector{Statement}(), true, region.entry, 1) :
        InsertBeforeLoopHead(Vector{Statement}(), region.loop, true)
    push!(call_sites.actions, action_before_region)
    
    for fknob in call_sites.extra.function_knobs
        mknobs = call_sites.extra.fknob2mknobs[fknob]
        for mknob in mknobs
            create_add_mknob_to_fknob(action_before_region.new_stmts,
                                      mknob, fknob)
        end
    end
end

@doc """"
Add function knob to the expressions (library function calls) in the AST.
"""
function add_function_knobs_to_calls(
    call_sites :: CallSites
)
    assert(!call_sites.extra.ast_may_change)

    for fknob in call_sites.extra.function_knobs
        expr = call_sites.extra.fknob2expr[fknob]
        push!(expr.args, fknob)
        
        # It seems that Julia looks up a dictionary's key by value, not by
        # address (It does not have the concept of pointer anyway), and thus
        # expr is considered changed with the new fknob added. So update the
        # mappings between the expr and the fknob:
        call_sites.extra.expr2fknob[expr]  = fknob
        call_sites.extra.fknob2expr[fknob] = expr
    end
end

@doc """
Create statements to propagate matrix info from one matrix knob to another if
the expression is an assignment statement.
"""
function propagate_matrix_info(
    expr        :: Any,
    call_sites :: CallSites
)
    assert(typeof(expr) == Expr)

    bb           = call_sites.extra.bb
    stmt_idx     = call_sites.extra.stmt_idx
    matrix_knobs = call_sites.extra.matrix_knobs
    matrix2mknob = call_sites.extra.matrix2mknob
    if expr.head == :(=) && length(expr.args) == 2
        # TODO: enable propagation for other kinds of assignment like >>=
        if haskey(matrix2mknob, expr.args[1]) && haskey(matrix2mknob, expr.args[2])
            mknob1 = matrix2mknob[expr.args[1]]
            mknob2 = matrix2mknob[expr.args[2]]
            create_propagate_matrix_info(mknob1, mknob2, bb, stmt_idx, call_sites)
        end
    end
    return nothing
end

@doc """ A handler to visit_expressions(). """
propagate_matrix_info(ast, call_sites :: CallSites, top_level_number, is_top_level, read) =
    propagate_matrix_info(ast, call_sites)

@doc """ 
Discover context-sensitive function calls in the loop region. Insert matrix and 
function-specific context info (mknobs and fknobs).
"""
function context_sensitive_transformation(
    actions         :: Vector{Action},
    region          :: Region,
    bb2depth        :: Dict{BasicBlock, Int}, 
    lambda          :: LambdaInfo,
    symbol_info     :: Sym2TypeMap,
    liveness        :: Liveness,
    cfg             :: CFG
)
    # We need to walk the AST recursively twice: first, match/replace and gather
    # context info. We cannot add context info into the AST yet, because
    # replacement happens recursively: a subtree might be replaced, and then
    # its parent tree is replaced subsequently, and during which the subtree is
    # removed. If we add context info into the subtree, some of its effect (the
    # actions to create them) remains, even though the subtree is gone. 
    # Another reason that we should not add context info in this step is that if
    # we do it, a subtree contains fknob, and then matching its parent tree needs
    # a pattern containing a fknob. That is rather problematic.
    # Therefore, we'd better walk the AST twice. After the AST is stablized after
    # the first walk, the second walk can add the context info.
    dprintln(1, 0, "\nContext-sensitive transformation:")
    recursive      = true
    ast_may_change = true
    call_sites = CallSites(Set{CallSite}(), region, lambda, symbol_info, liveness,
                           CS_transformation_patterns,
                           actions, ContextExtra(region.symbol_property, bb2depth, ast_may_change))
    
#    global trace_call_replacement = true
    
    # All patterns are non-splittable.
    call_sites.extra.non_splitting_patterns = copy(CS_transformation_patterns)
    visit_expressions(region, cfg, call_sites, recursive, match_replace_an_expr_pattern, match_replace_an_two_statements_pattern)

#    global trace_call_replacement = false

    # The second walk: gather the knobs of the nodes in the AST. The knobs
    # were generated in the first walk. Not all those generated are gathered,
    # since some ast nodes have been replaced, and thus won't be visited this
    # time.
    call_sites.extra.ast_may_change = false
    visit_expressions(region, cfg, call_sites, recursive, gather_knobs)

    # Generate/delete the knobs at the beginning/exit of the loop
    generate_and_delete_knobs(region, call_sites)

    # Now at the beginning of the loop, associate function knobs with their 
    # matrix knobs.
    add_matrix_knobs_to_function_knobs(call_sites)

    # Propagate matrix info from one matrix knob to another for assignment statements
    recursive  = false
    visit_expressions(region, cfg, call_sites, recursive, propagate_matrix_info)

    # Add fknobs to function calls in the AST. This changes AST, and thus should 
    # be put at the last, because otherwise, the maps in the extra including
    # matrix2mknobs and expr2fknob whose keys are AST nodes, might become
    # out of date.
    add_function_knobs_to_calls(call_sites)

    # Reordering is now part of context-sensitive optimization.
    if reorder_enabled
        reordering(actions, region, symbol_info, liveness, cfg, call_sites)
    end
    
    #global trace_call_replacement = false
end

@doc """ 
Discover context-sensitive function calls in the loop regions. Insert matrix and 
function-specific context info (mknobs and fknobs) into actions.
"""
function AST_context_sensitive_transformation(
    actions     :: Vector{Action},
    func_region :: FunctionRegion,
    regions     :: Vector{LoopRegion},
    bb2depth    :: Dict{BasicBlock, Int},
    lambda      :: LambdaInfo,
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG
)
    # Find properties for the function and loop regions.
    func_region.symbol_property = find_properties_of_matrices(func_region, lambda, symbol_info, liveness, cfg)
    for region in regions
        region.symbol_property = find_properties_of_matrices(region, lambda, symbol_info, liveness, cfg)
    end

    if context_sensitive_opt_for_func
        # Context opt for the function
        context_sensitive_transformation(actions, func_region, bb2depth, lambda, symbol_info, liveness, cfg)
    else
        # Context opt for each loop region
       for region in regions
            context_sensitive_transformation(actions, region, bb2depth, lambda, symbol_info, liveness, cfg)
       end
    end
    
    dprintln(1, 0, "\nContext-sensitive actions to take:")
    dprintln(1, 1, "", actions)
    
    actions
end
