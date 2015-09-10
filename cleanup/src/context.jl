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
    # TODO: replace the code to be based on call_sites.matrix_properties.

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
    # So far, we handle only a loop region
    assert(typeof(call_sites.region) == LoopRegion)

    # Create a function-specific knob at each call site of the context-specific
    # functions. Create matrix-specific knobs for the matrices inputs.
    region               = call_sites.region
    L                    = region.loop
    symbol_info          = call_sites.symbol_info
    site                 = CallSite(ast)
    args                 = ast.args
    action_before_region = InsertBeforeLoopHead(Vector{Statement}(), L, true)
    fknob                = create_new_function_knob(action_before_region.new_stmts,
                                                    ast, fknob_creator)
    push!(call_sites.actions, action_before_region)
    push!(call_sites.sites, site)
    for arg in matrices_to_track
        if arg == :result
            M = ast
        else
            new_arg = replacement_arg(arg, args, nothing, symbol_info)
            assert(typeof(new_arg) == Symbol || typeof(new_arg) == SymbolNode ||
                   typeof(new_arg) == GenSym)
            M = ((typeof(new_arg) == SymbolNode) ? new_arg.name : new_arg)
        end
        if !haskey(call_sites.matrix_knobs, M)
            # Create statements that will create and intialize a knob for
            # the matrix before the loop region
            mknob = create_new_matrix_knob(action_before_region.new_stmts, M,
                                           call_sites)
            call_sites.matrix_knobs[M] = mknob
        end
        create_add_mknob_to_fknob(action_before_region.new_stmts,
                                  call_sites.matrix_knobs[M], fknob)
    end

    # Create statements that will delete the function knob at region exits
    for exit in region.exits
        action  = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
        create_delete_function_knob(action.new_stmts, fknob_deletor, fknob)
        push!(call_sites.actions, action)
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
Post-processing function of CS_ADAT_pattern
"""
function CS_ADAT_post_replacement(
    ast               :: Expr,
    call_sites        :: CallSites,
    fknob_creator     :: String,
    fknob_deletor     :: String,
    matrices_to_track :: Tuple
)
    # We need to replace A * D * A' into ADB(A', D, A, fknob). At this 
    # moment, it is in the middle form of ADB(A, D, A). That means:
    # (1) Make a symbol AT. Insert AT = A' before the loop (This is to hoist A'
    #     out of loop)
    # (2) Replace ADB(A, D, A) as CSR_ADB(AT, D, A).
    # fknob is not added for now, which will be added automatically later.
    action = InsertBeforeLoopHead(Vector{Statement}(), call_sites.region.loop, true)
    push!(call_sites.actions, action)

    A  = ast.args[2]
    assert(typeof(A) == SymbolNode)
    AT = Symbol(string("__", string(A.name), "T__"))
    stmt = Statement(-1, Expr(:(=), AT, Expr(:call, GlobalRef(Main, :ctranspose), A)))
    push!(action.new_stmts, stmt)
    
    ast.args[2] = AT

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
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:ADB)),
     :arg2, :arg3, :arg2),
    CS_ADAT_post_replacement,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:result, :arg2, :arg3, :arg4)
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
    (:result, :arg2)
)

const CS_cholsolve_pattern = ExprPattern(
    "CS_cholsolve_pattern",
    (:call, GlobalRef(Main, :\), Base.SparseMatrix.CHOLMOD.Factor{Float64}, Any),
    (:NO_SUB_PATTERNS,),
    do_nothing,
    (:call, TypedExprNode(Function, :call, TopNode(:getfield), :SparseAccelerator, QuoteNode(:cholmod_factor_inverse_divide)),
     :arg2, :arg3),
    gather_context_sensitive_info,
    "NewFunctionKnob",
    "DeleteFunctionKnob",
    (:arg2,)
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
    CS_ADAT_pattern,
    CS_cholfact_int32_pattern,
    CS_cholsolve_pattern,
    CS_fwdTriSolve!_pattern,
    CS_bwdTriSolve!_pattern
]

@doc """
Create statements that will create a matrix knob for matrix M.
"""
function create_new_matrix_knob(
    new_stmts  :: Vector{Statement},
    M          :: Symexpr,
    call_sites :: CallSites
)
    if !haskey(call_sites.matrix_properties, M)
        # This is a hack! Just to unblock us temporarily.
        # TODO: Once Todd's structure analysis is done, this hack MUST be removed!
        call_sites.matrix_properties[M] = MatrixProperties(false, true, false, false, false, false)
        # TODO: Enable this once the hack is removed.
        #throw(MatrixPropertiesUnavail(M))
    end
    properties = call_sites.matrix_properties[M]
    mknob      = gensym(string("mknob", typeof(M) == Expr ? "Expr" : string(M)))

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
        # This is a hack, just to unblock us temporarily.
        # TODO: once Todd's structure analysis is done, this MUST be removed!
        if !properties.constant_structured
            call_sites.matrix_properties[M].constant_structured = true
            properties.constant_structured = true
        end

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
    mknob
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
                     Expr(:call, GlobalRef(SparseAccelerator, :PropagateMatrixInfo), 
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
    ast           :: Expr,
    fknob_creator :: String,
)
    assert(fknob_creator != "")
    
    fknob = gensym("fknob")
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
    ast.args = [ast.args; fknob]

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
    L                 = region.loop
    blocks            = cfg.basic_blocks
    matrix_properties = find_properties_of_matrices(region, liveness, cfg)
    matrix_knobs      = Dict{Symexpr, Symbol}()
    call_sites        = CallSites(Set{CallSite}(), region, symbol_info, 
                            matrix_properties, CS_transformation_patterns,
                            actions, Dict{Symexpr, Symbol}())
    for bb_idx in L.members
        bb         = blocks[bb_idx]
        statements = bb.statements
        for stmt_idx in 1 : length(statements)
            stmt = statements[stmt_idx]
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
            end

            # Create a function-specific knob at each call site of the context-specific
            # functions. Create matrix-specific knobs for the matrices inputs.
            # Create statements that will create and intialize a knob for
            # the matrix before the loop region, and statements that will delete
            # all the function knobs at each exit of the region
            CompilerTools.AstWalker.AstWalk(expr, match_replace_an_expr_pattern, call_sites)
        end
    end

    # We have modified call_sites.matrix_knobs. Now pass results back
    matrix_knobs = call_sites.matrix_knobs

    # Propagate matrix info from one to another for assignment statements
    for bb_idx in L.members
        bb         = blocks[bb_idx]
        statements = bb.statements
        for stmt_idx in 1 : length(statements)
            stmt = statements[stmt_idx]
            expr = stmt.expr
            if typeof(expr) != Expr
                continue
             end

            if expr.head == :(=) && length(expr.args) == 2
                # TODO: enable propagation for other kinds of assignment like >>=
                if haskey(matrix_knobs, expr.args[1]) && 
                   haskey(matrix_knobs, expr.args[2])
                   mknob1 = matrix_knobs[expr.args[1]]
                   mknob2 = matrix_knobs[expr.args[2]]
                   create_propagate_matrix_info(mknob1, mknob2, bb, stmt_idx, call_sites)
                end
            end
        end
    end

    # Create statements that will delete the matrix knobs at each exit of the region
    for mknob in values(matrix_knobs)
        for exit in region.exits
            action = InsertOnEdge(Vector{Statement}(), exit.from_bb, exit.to_bb)
            create_delete_matrix_knob(action.new_stmts, mknob)
            push!(call_sites.actions, action)
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