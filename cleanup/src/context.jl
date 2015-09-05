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
    matrices_to_track :: Tuple
)
    symbol_info = call_sites.symbol_info
    site        = CallSite(ast, Vector(), Vector(), fknob_creator, fknob_deletor)
    args        = ast.args
    for arg in matrices_to_track
        new_arg = replacement_arg(arg, args, nothing, symbol_info)
        assert(typeof(new_arg) == Symbol || typeof(new_arg) == SymbolNode ||
               typeof(new_arg) == GenSym)
        push!(site.matrices, typeof(new_arg) == SymbolNode ? new_arg.name : new_arg)
    end
    push!(call_sites.sites, site)
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

    # TODO: enable this once liveness def() bug is fixed
#    if !in(A, call_sites.constants)
#    println("A not const: ", A)
#    println("constants = ", call_sites.constants)
#        return false
#    end
    
    # TODO: check D is diagonal after structure analysis is done.
    return true
end

@doc """ 
Post-processing function of CS_ADAT_assign_pattern
"""
function CS_ADAT_assign_post_replacement(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple
)
    # We need to replace C = A * D * A' into ADB(C, A', D, A, fknob). At this 
    # moment, it is in the middle form of ADB(C, A, D, A). That means:
    # (1) Make a symbol AT. Insert AT = A' before the loop (This is to hoist A'
    #     out of loop)
    # (2) Replace ADB(C, A, D, A) as CSR_ADB(C, AT, D, A).
    # fknob is not added for now, which will be added automatically later.
    action = InsertBeforeLoopHead(Vector{Statement}(), call_sites.region.loop, true)
    push!(call_sites.actions, action)

    A  = ast.args[3]
    assert(typeof(A) == SymbolNode)
    AT = Symbol(string("__", string(A.name), "T__"))
    stmt = Statement(-1, Expr(:(=), AT, Expr(:call, GlobalRef(Main, :ctranspose), A)))
    push!(action.new_stmts, stmt)
    
    ast.args[3] = AT

    return gather_context_sensitive_info(ast, call_sites, fknob_creator, fknob_deletor, matrices_to_track)
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
    ()
)

const CS_ADAT_pattern = ExprPattern(
    "CS_ADAT_pattern",
    (:call, GlobalRef(Main, :*), SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC),
    (nothing, nothing, nothing, nothing, CS_ADAT_AT_pattern),
    CS_ADAT_check,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const CS_ADAT_assign_pattern = ExprPattern(
    "CS_ADAT_assign_pattern",
    (:(=), SparseMatrixCSC, SparseMatrixCSC),
    (nothing, nothing, CS_ADAT_pattern),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:ADB)),
      :arg1, :aarg22, :aarg23, :aarg22),
    CS_ADAT_assign_post_replacement,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2, :arg3, :arg4, :arg5)
)

const CS_cholfact_int32_pattern = ExprPattern(
    "CS_cholfact_int32_pattern",
    (:call, GlobalRef(Main, :cholfact_int32), SparseMatrixCSC{Float64, Int32}),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const CS_cholfact_int32_assign_pattern = ExprPattern(
    "CS_cholfact_int32_assign_pattern",
    (:(=), Base.SparseMatrix.CHOLMOD.Factor{Float64}, Base.SparseMatrix.CHOLMOD.Factor{Float64}),
    (nothing, nothing, CS_cholfact_int32_pattern),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:cholfact_int32)),
     :arg1, :aarg22),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2, :arg3)
)

const CS_cholsolve_pattern = ExprPattern(
    "CS_cholsolve_pattern",
    (:call, GlobalRef(Main, :\), Base.SparseMatrix.CHOLMOD.Factor{Float64}, Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:NO_CHANGE, ),
    do_nothing,
    "",
    "",
    ()
)

const CS_cholsolve_assign_pattern = ExprPattern(
    "CS_cholsolve_assign_pattern",
    (:(=), Any, Any),
    (nothing, nothing, CS_cholsolve_pattern),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:cholmod_factor_inverse_divide)),
     :arg1, :aarg22, :aarg23),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg3,)
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
    (:arg2,)
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
    (:arg2,)
)

@doc """" Patterns that will actually transform the code. """
CS_transformation_patterns = [
    CS_ADAT_assign_pattern,
    CS_cholfact_int32_assign_pattern,
    CS_cholsolve_assign_pattern,
    CS_fwdTriSolve!_pattern,
    CS_bwdTriSolve!_pattern
]

@doc """
Create statements that will create a matrix knob for matrix M.
"""
function create_new_matrix_knob(
    new_stmts :: Vector{Statement},
    M         :: Union{Sym, Expr}
)
    mknob = gensym(string("mknob", string(M)))
    
    #TODO: set the right values here
    constant_valued        = false
    constant_structured    = false
    is_symmetric           = false
    is_structure_symmetric = false
    is_structure_only      = false

    new_stmt = Expr(:(=), mknob,
                Expr(:call, GlobalRef(SparseAccelerator, :new_matrix_knob), 
                     M, constant_valued, constant_structured, is_symmetric,
                     is_structure_symmetric, is_structure_only))
    push!(new_stmts, Statement(0, new_stmt))
    
    mknob
end

@doc """
Create statements that will increment a matrix knob's version.
"""
function create_increment_matrix_version(
    new_stmts :: Vector{Statement},
    mknob     :: Symbol
)
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :increment_matrix_version), mknob)
    push!(new_stmts, Statement(0, new_stmt))
end

@doc """
Create statements that will create a function knob for the call site, and add
the function knob to the call as a parameter.
"""
function create_new_function_knob(
    new_stmts :: Vector{Statement},
    call_site :: CallSite
)
    assert(call_site.fknob_creator != "")
    
    fknob = gensym("fknob")
    if call_site.fknob_creator == "NewFunctionKnob"
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
                              call_site.fknob_creator
                        )
                   )
    end
    push!(new_stmts, Statement(0, new_stmt))
    call_site.ast.args = [call_site.ast.args; fknob]

    fknob
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
    fknob_deletor :: String,
    fknob         :: Symbol
)
    assert(call_site.fknob_deletor != "")

    if call_site.fknob_deletor == "DeleteFunctionKnob"
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
    L           = region.loop
    blocks      = cfg.basic_blocks
    constants   = find_constant_values(region, liveness, cfg)
    call_sites  = CallSites(Set{CallSite}(), region, symbol_info, constants, CS_transformation_patterns, actions)
    var_defs    = Dict{Sym, Set{Tuple{BasicBlock, StatementIndex}}}() # Map from a variable to a set of statements defining it
    for bb_idx in L.members
        bb         = blocks[bb_idx]
        statements = bb.statements
        for stmt_idx in 1 : length(statements)
            stmt = statements[stmt_idx]
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)

            stmt_def = LivenessAnalysis.def(stmt, liveness)
            for d in stmt_def
                if type_of_ast_node(d, symbol_info) <: AbstractSparseMatrix
                    if !haskey(var_defs, d)
                        var_defs[d] = Set{Tuple{BasicBlock, StatementIndex}}()
                    end
                    push!(var_defs[d], (bb, stmt_idx))
                end
            end
        end
    end

    # Create a function-specific knob at each call site of the context-specific
    # functions. Create matrix-specific knobs for the matrices inputs.
    # First, create matrix knobs, as they will be needed for creating the function
    # knobs.
    matrix_knobs         = Dict{Any, Symbol}()
    action_before_region = InsertBeforeLoopHead(Vector{Statement}(), L, true)
    push!(actions, action_before_region)
    for call_site in call_sites.sites
        for M in call_site.matrices
            if !haskey(matrix_knobs, M)
                # Create statements that will create and intialize a knob for
                # the matrix before the loop region
                mknob = create_new_matrix_knob(action_before_region.new_stmts, M)
                matrix_knobs[M] = mknob
            end

            if in(M, call_site.matrices_to_track_values)
                # Create statements that will update the knob before every
                # statement that defines the matrix
                if haskey(var_defs, M)
                    for (bb, stmt_idx) in var_defs[M]
                        action = InsertBeforeStatement(Vector{Statement}(), bb, stmt_idx)
                        push!(actions, action)
                        mknob = matrix_knobs[M]
                        create_increment_matrix_version(action.new_stmts, mknob)
                    end
                end
            end
        end
    end

    function_knobs = Set()
    for call_site in call_sites.sites
        fknob = create_new_function_knob(action_before_region.new_stmts, call_site)
        push!(function_knobs, (fknob, call_site.fknob_deletor))
        for M in call_site.matrices
            create_add_mknob_to_fknob(action_before_region.new_stmts, matrix_knobs[M], fknob)
        end
    end
    
    # Create statemetns that will delete all the knobs at each exit of the region
    for exit in region.exits
        action  = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
        push!(actions, action)
        for (fknob, fknob_deletor) in function_knobs
            create_delete_function_knob(action.new_stmts, fknob_deletor, fknob)
        end
        for mknob in values(matrix_knobs)
            create_delete_matrix_knob(action.new_stmts, mknob)
        end
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
    # Discover structures of matrices in the whole function
    structure_discovery(symbol_info, cfg)

    # Discover context info for each loop region
    for region in regions
        context_info_discovery(actions, region, func_ast, symbol_info, liveness, cfg)
    end

    dprintln(1, 0, "\nContext-sensitive actions to take:")
    dprintln(1, 1, "", actions)
    
    actions
end