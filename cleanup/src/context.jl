@doc """
The CallSites' extra field for context info discovery.
"""
type ContextExtra
    # Environment before the AST walking
    matrix_properties        :: Symexpr2PropertiesMap

    # Flag if ast may be under change. We visite the AST twice. It is true
    # the first time (allowing pattern replacement; knobs are generated as well),
    # and false the second time (allowing gathering the knobs of the AST nodes).
    ast_may_change           :: Bool

    # Infomration gathered when AST is walked the first time
    reordering_decider       :: Any    # It should be expr, but then cannot be initialized with nothing
    reordering_decider_power :: Int    # 0 is no power: not to make reordering decision
    reordering_FAR           :: Vector # First Arrays Reordered by the decider
    expr2fknob               :: Dict{Expr, Symbol}
    fknob2expr               :: Dict{Symbol, Expr}
    fknob2mknobs             :: Dict{Symbol, Vector} # Vector[matrix knobs]
    fknob2ctor_dtor          :: Dict{Symbol, Tuple}  # Tuple(fknob_creator, fknob_deletor)
    matrix2mknob             :: Dict{Symexpr, Symbol}
    mknob2matrix             :: Dict{Symbol, Symexpr}

    # The knobs gathered when AST is walked the second time.
    matrix_knobs             :: Set{Symbol}
    function_knobs           :: Set{Symbol}
 
    # Scratch variables.
    bb                       :: Any              # BasicBlock
    stmt_idx                 :: StatementIndex
    
    ContextExtra(_matrix_properties, _ast_may_change) = new(
            _matrix_properties, _ast_may_change, nothing, 0, [],
            Dict{Expr, Symbol}(), Dict{Symbol, Expr}(), Dict{Symbol, Vector}(),
            Dict{Symbol, Tuple}(), Dict{Symexpr, Symbol}(),
            Dict{Symbol, Symexpr}(),
            Set{Symbol}(), Set{Symbol}(), nothing, 0)
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
    fknob_creator :: String,
    fknob_deletor :: String
)
    # TODO: replace the code to be based on call_sites.extra.matrix_properties.

    L_or_U = ast.args[2]
    structure = get_structure_proxy(L_or_U) 
    if structure == nothing
        return false
    end
    
    # Check it is lower or upper part of a symmetric matrix.
    if structure.lower || structure.upper
        # So far, symmetricity is not checked: that might need user annotation.
        # TODO: Check symmetricity once available
        return true
    end
    return false
end

@doc """ 
Post-processing function of a pattern: Gather context sensitive info of the call
site at the given AST.
"""
function gather_context_sensitive_info(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple,
    reordering_power  :: Int,
    reordering_FAR    :: Tuple
)
    assert(call_sites.extra.ast_may_change)

    # Create a function-specific knob. 
    fknob                                   = gensym("fknob")
    call_sites.extra.expr2fknob[ast]        = fknob
    call_sites.extra.fknob2expr[fknob]      = ast
    call_sites.extra.fknob2ctor_dtor[fknob] = (fknob_creator, fknob_deletor)
    call_sites.extra.fknob2mknobs[fknob]    = []
    
    # Create matrix-specific knobs for the matrices inputs.
    symbol_info = call_sites.symbol_info
    args        = ast.args
    for arg in matrices_to_track
        if arg == :result
            M = ast
        else
            new_arg = replacement_arg(arg, args, nothing, symbol_info)
            assert(typeof(new_arg) == Symbol || typeof(new_arg) == SymbolNode ||
                   typeof(new_arg) == GenSym)
            M = ((typeof(new_arg) == SymbolNode) ? new_arg.name : new_arg)
        end
        if !haskey(call_sites.extra.matrix2mknob, M)
            mknob = gensym(string("mknob", typeof(M) == Expr ? "Expr" : string(M)))
            call_sites.extra.matrix2mknob[M]     = mknob
            call_sites.extra.mknob2matrix[mknob] = M
        end
        push!(call_sites.extra.fknob2mknobs[fknob], call_sites.extra.matrix2mknob[M])
    end

    if call_sites.extra.reordering_decider_power < reordering_power
        call_sites.extra.reordering_decider       = ast
        call_sites.extra.reordering_decider_power = reordering_power
        call_sites.extra.reordering_FAR           = []
        for arg in reordering_FAR
            new_arg = replacement_arg(arg, args, nothing, symbol_info)
            assert(typeof(new_arg) == Symbol || typeof(new_arg) == SymbolNode ||
                   typeof(new_arg) == GenSym)
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
    fknob_creator :: String,
    fknob_deletor :: String
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
    fknob_creator     :: String,
    fknob_deletor     :: String,
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
    action = InsertBeforeLoopHead(Vector{Statement}(), call_sites.region.loop, true)
    push!(call_sites.actions, action)

    A = ast.args[4]
    assert(typeof(A) == SymbolNode)
    properties = call_sites.extra.matrix_properties[A.name]
    AT = Symbol(string("__", string(A.name), "T__"))
    stmt = Statement(-1, Expr(:(=), AT, Expr(:call, GlobalRef(Main, :ctranspose), A)))
    push!(action.new_stmts, stmt)

    ast.args[4] = SymbolNode(AT, A.typ)

    # AT has the same properties as A. But AT has only a single static definition.
    call_sites.extra.matrix_properties[AT]               = properties
    call_sites.extra.matrix_properties[AT].is_single_def = true

    # The result is a symmetric matrix.
    # ISSUE: in analysis like constant value/structure analysis, should we 
    # consider each AST node instead of only symbols? If that is the case, we
    # do not have to make properties here.
    call_sites.extra.matrix_properties[ast]                        = MatrixProperties()
    call_sites.extra.matrix_properties[ast].constant_valued        = false #TODO: should determined by if A and D are constant valued
    call_sites.extra.matrix_properties[ast].constant_structured    = true
    call_sites.extra.matrix_properties[ast].is_symmetric           = true
    call_sites.extra.matrix_properties[ast].is_structure_symmetric = true
    call_sites.extra.matrix_properties[ast].is_structure_only      = false
    call_sites.extra.matrix_properties[ast].is_single_def          = true

    return gather_context_sensitive_info(ast, call_sites, fknob_creator, 
                                         fknob_deletor, matrices_to_track,
                                         reordering_power, reordering_FAR)
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
        # SparseAccelerator.ADB(A', D, A, fknob_ADB)
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
        B = SparseAccelerator.ADB(A', D, A, fknob_ADB)
        ...
        # SparseAccelerator.cholfact_int32(B, fknobcholfact_int32)
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
        
        R = SparseAccelerator.cholfact_int32(B, fknob_cholfact_int32)
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

const CS_ADAT_AT_pattern = ExprPattern(
    "CS_ADAT_AT_pattern",
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
    (nothing, nothing, nothing, nothing, CS_ADAT_AT_pattern),
    CS_ADAT_check,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:ADB)),
     :arg2, :arg3, :arg2),
    CS_ADAT_post_replacement,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:result, :arg2, :arg3, :arg4),
    0,
    ()
)

const CS_cholfact_int32_pattern = ExprPattern(
    "CS_cholfact_int32_pattern",
    (:call, GlobalRef(Main, :cholfact_int32), SparseMatrixCSC{Float64, Int32}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:cholfact_int32)),
     :arg2),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2,),
    0,
    ()
)

const CS_cholsolve_pattern = ExprPattern(
    "CS_cholsolve_pattern",
    (:call, GlobalRef(Main, :\), Base.SparseMatrix.CHOLMOD.Factor{Float64}, Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:cholfact_inverse_divide)),
     :arg2, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2,),
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
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:fwdTriSolve!)),
      :arg2, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2,),
    100,
    (:arg2, :arg3)
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
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:bwdTriSolve!)),
      :arg2, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2,),
    90,
    (:arg2, :arg3)
)

@doc """
The following patterns are for SpMV. The difference from the SpMV patterns in
call-replacement.jl is that here we will add mknobs and fknob to SpMV.
"""

@doc """ A * x => SpMV(A, x) """
const CS_SpMV_pattern1 = ExprPattern(
    "CS_SpMV_pattern1",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2,),
    40,
    (:arg2, :arg3)
)

@doc """ a * A * x + g => SpMV(a, A, x, 0, x, g) """
const CS_SpMV_pattern2 = ExprPattern(
    "CS_SpMV_pattern2",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, number_times_matrix_vector_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg3,),
    40,
    (:arg3, :arg4)
)

@doc """ a * A * x => SpMV(a, A, x) """
const CS_SpMV_pattern3 = ExprPattern(
    "CS_SpMV_pattern3",
    (:call, GlobalRef(Main, :*), Number, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :arg2, :arg3, :arg4),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg3,),
    40,
    (:arg3, :arg4)
)

@doc """ SpMV(a, A, x) + g => SpMV(a, A, x, 0, x, g) """
const CS_SpMV_pattern4 = ExprPattern(
    "CS_SpMV_pattern4",
    (:call, GlobalRef(Main, :+), Vector, Number),
    (nothing, nothing, SpMV_3_parameters_pattern, nothing),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV)),
     :aarg22, :aarg23, :aarg24, 0, :aarg24, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg3,),
    40,
    (:arg3,)
)

@doc """ A_mul_B!(y, A, x) = SpMV!(y, A, x) """
const CS_SpMV!_pattern1 = ExprPattern(
    "CS_SpMV!_pattern1",
    (:call, GlobalRef(Main, :A_mul_B!), Vector, SparseMatrixCSC, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg2, :arg3, :arg4),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg3,),
    50,
    (:arg3, :arg4, :arg2) # Put seed (A) as the first
)

@doc """ z = SpMV(a, A, x, b, y, g), z is x or y => SpMV!(z, a, A, x, b, y, g) """
const CS_SpMV!_pattern2 = ExprPattern(
    "CS_SpMV!_pattern2",
    (:(=), Vector, Vector),
    (nothing, nothing, SpMV_6_parameters_pattern),
    LHS_in_RHS,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg1, :aarg22, :aarg23, :aarg24, :aarg25, :aarg26, :aarg27),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg4,),
    80,
    (:arg4, :arg2, :arg5, :arg7) # Put seed (A) at the first
)

@doc """ x = a * A * x => SpMV!(x, a, A, x, 0, x, 0) """
const CS_SpMV!_pattern3 = ExprPattern(
    "CS_SpMV!_pattern3",
    (:(=), Vector, Vector),
    (nothing, nothing, SpMV_3_parameters_pattern),
    LHS_in_RHS,
    (TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg1, :aarg22, :aarg23, :aarg24, 0.0, :arg1, 0.0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg4,),
    60,
    (:arg4, :arg2, :arg5) # Put seed (A) at the first
)

@doc """ A_mul_B!(a, A, x, b, y) => SpMV!(y, a, A, x, b, y , 0) """
const CS_SpMV!_pattern4 = ExprPattern(
    "CS_SpMV!_pattern4",
    (:call, GlobalRef(Main, :A_mul_B!), Number, SparseMatrixCSC, Vector, Number, Vector),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:SpMV!)),
     :arg6, :arg2, :arg3, :arg4, :arg5, :arg6, 0.0),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg4,),
    70,
    (:arg4, :arg2, :arg5, :arg7) # Put seed (A) at the first
)

@doc """" Patterns that will actually transform the code. """
CS_transformation_patterns = [
    CS_ADAT_pattern,
    CS_cholfact_int32_pattern,
    CS_cholsolve_pattern,
    CS_fwdTriSolve!_pattern,
    CS_bwdTriSolve!_pattern,
    CS_SpMV_pattern1,
    CS_SpMV_pattern2,
    CS_SpMV_pattern3,
    CS_SpMV_pattern4,
    CS_SpMV!_pattern1,
    CS_SpMV!_pattern2,
    CS_SpMV!_pattern3,
    CS_SpMV!_pattern4
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
    if properties.constant_valued
        new_stmt = Expr(:(=), mknob,
                    Expr(:call, GlobalRef(SparseAccelerator, :new_matrix_knob), 
                         M,
                         properties.constant_valued, 
                         properties.constant_structured, 
                         properties.is_symmetric,
                         properties.is_structure_symmetric,
                         properties.is_structure_only,
                         properties.is_single_def))
    else
        assert(properties.constant_structured) 
        new_stmt = Expr(:(=), mknob,
                    Expr(:call, GlobalRef(SparseAccelerator, :new_matrix_knob), 
                         properties.constant_valued, 
                         properties.constant_structured, 
                         properties.is_symmetric,
                         properties.is_structure_symmetric,
                         properties.is_structure_only,
                         properties.is_single_def))
    end
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
    stmt = Statement(-1, 
                     Expr(:call, GlobalRef(SparseAccelerator, :propagate_matrix_info), 
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
    fknob_creator :: String,
)
    assert(fknob_creator != "")

    if fknob_creator == "NewFunctionKnob"
        # Create a function knob. It has no private info specific to that function.
        # Call the parameterless version of new_function_knob for cleaner code.
        new_stmt = Expr(:(=), fknob, 
                        Expr(:call, GlobalRef(SparseAccelerator, :new_function_knob))
                   )
    else
        # So far, we do not create any special knob. So this case is never reached.
        assert(false)
        new_stmt = Expr(:(=), fknob, 
                        Expr(:call, GlobalRef(SparseAccelerator, :new_function_knob),
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
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :add_mknob_to_fknob), mknob, fknob)
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create statements that will delete the function knob.
"""
function create_delete_function_knob(
    new_stmts     :: Vector{Statement},
    fknob         :: Symbol,
    fknob_deletor :: String
)
    assert(fknob_deletor != "")

    if fknob_deletor == "DeleteFunctionKnob"
        # Delete a function knob, which has no private info specific to that function.
        # Call the 1-parameter version of delete_function_knob for cleaner code.
        new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :delete_function_knob), 
                         fknob)
    else
        # So far, we have not created any special knob. So this case is never reached.
        assert(false)
        new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :delete_function_knob),
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
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :delete_matrix_knob), mknob)
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Visit each expression the loop region. If recursive is true, visit each AST node
of the expression and handle each with the handler. Otherwise, handle the
expression with the handler.
"""
function vist_expressions(
    region      :: LoopRegion,
    cfg         :: CFG,
    call_sites  :: CallSites,
    recursive   :: Bool,
    handler     :: Function
)
    L      = region.loop
    blocks = cfg.basic_blocks
    for bb_idx in L.members
        bb                  = blocks[bb_idx]
        statements          = bb.statements
        call_sites.extra.bb = bb
        for stmt_idx in 1 : length(statements)
            call_sites.extra.stmt_idx = stmt_idx
            stmt                      = statements[stmt_idx]
            expr                      = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            if recursive
                CompilerTools.AstWalker.AstWalk(expr, handler, call_sites)
            else
                handler(expr, call_sites)
            end
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

@doc """
Create statements that will create and intialize a knob for a matrix/function
before the loop, and statements that will delete all the matrix/function
knobs at each exit of the loop.
"""
function generate_and_delete_knobs(
    region     :: LoopRegion,
    call_sites :: CallSites
)
    # So far, we handle only a loop region
    assert(typeof(call_sites.region) == LoopRegion)

    region         = call_sites.region
    L              = region.loop
    matrix_knobs   = call_sites.extra.matrix_knobs
    function_knobs = call_sites.extra.function_knobs
    symbol_info    = call_sites.symbol_info

    action_before_region = InsertBeforeLoopHead(Vector{Statement}(), L, true)
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

    # Create statements that will delete the knobs at region exits
    for exit in region.exits
        action  = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
        push!(call_sites.actions, action)

        for mknob in matrix_knobs
            create_delete_matrix_knob(action.new_stmts, mknob)
        end

        for fknob in function_knobs
            fknob_creator, fknob_deletor = call_sites.extra.fknob2ctor_dtor[fknob]
            create_delete_function_knob(action.new_stmts, fknob, fknob_deletor)
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

   # So far, we handle only a loop region
    assert(typeof(call_sites.region) == LoopRegion)

    L                    = call_sites.region.loop
    action_before_region = InsertBeforeLoopHead(Vector{Statement}(), L, true)
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
function-specific context info (mknobs and fknobs) into actions.
"""
function context_info_discovery(
    actions     :: Vector{Action},
    region      :: LoopRegion,
    func_ast    :: Expr,
    symbol_info :: Sym2TypeMap,
    liveness    :: Liveness,
    cfg         :: CFG
)
    matrix_properties = find_properties_of_matrices(region, symbol_info, liveness, cfg)

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
    recursive      = true
    ast_may_change = true
    call_sites = CallSites(Set{CallSite}(), region, symbol_info,
                           CS_transformation_patterns,
                           actions, ContextExtra(matrix_properties, ast_may_change))
    vist_expressions(region, cfg, call_sites, recursive, match_replace_an_expr_pattern)

    # The second walk: gather the knobs of the nodes in the AST. The knobs
    # were generated in the first walk. Not all those generated are gathered,
    # since some ast nodes have been replaced, and thus won't be visited this
    # time.
    call_sites.extra.ast_may_change = false
    vist_expressions(region, cfg, call_sites, recursive, gather_knobs)

    # Generate/delete the knobs at the beginning/exit of the loop
    generate_and_delete_knobs(region, call_sites)

    # Now at the beginning of the loop, associate function knobs with their 
    # matrix knobs.
    add_matrix_knobs_to_function_knobs(call_sites)

    # Add fknobs to function calls in the AST
    add_function_knobs_to_calls(call_sites)

    # Propagate matrix info from one matrix knob to another for assignment statements
    recursive  = false
    vist_expressions(region, cfg, call_sites, recursive, propagate_matrix_info)

    # Reordering is now part of context-sensitive optimization.
    if reorder_enabled
        reordering(actions, region, symbol_info, liveness, cfg, call_sites)
    end
end

@doc """ 
Discover context-sensitive function calls in the loop regions. Insert matrix and 
function-specific context info (mknobs and fknobs) into actions.
"""
function context_info_discovery(
    actions     :: Vector{Action},
    regions     :: Vector{LoopRegion},
    func_ast    :: Expr, 
    symbol_info :: Sym2TypeMap, 
    liveness    :: Liveness, 
    cfg         :: CFG
)
    # Discover context info for each loop region
    for region in regions
        context_info_discovery(actions, region, func_ast, symbol_info, liveness, cfg)
    end

    dprintln(1, 0, "\nContext-sensitive actions to take:")
    dprintln(1, 1, "", actions)
    
    actions
end
