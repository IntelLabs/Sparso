@doc """
For a function, describe the matrix arguments, and how to create and delete its
context info (fknob).
"""
immutable ContextSensitiveFunction
    module_name     :: String # Module of the function. 
    function_name   :: String # Name of the function
    argument_types  :: Tuple  # Tuple of the function arguments' types
    matrices        :: Set    # The matrix arguments.
    fknob_creator   :: Tuple  # {Symbol, String} The path to a fknob creator
    fknob_deletor   :: Tuple  # {Symbol, String} The path to a fknob deletor
end

# Below are the context sensitive functions we care about. For short, CS represents 
# Context Sensitive
const CS_fwdTriSolve! = ContextSensitiveFunction(
    "Base.SparseMatrix", 
    "fwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),
    (:NewForwardTriangularSolveKnob, libcsr),
    (:DeleteForwardTriangularSolveKnob, libcsr)
)

const CS_bwdTriSolve! = ContextSensitiveFunction(
    "Base.SparseMatrix", 
    "bwdTriSolve!",                              
    (AbstractSparseMatrix, Vector),
    Set(2),
    (:NewBackwardTriangularSolveKnob, libcsr),
    (:DeleteBackwardTriangularSolveKnob, libcsr)
)

context_sensitive_functions = [
    CS_fwdTriSolve!,
    CS_bwdTriSolve!
]

@doc """
Find a call to a context-sensitive function.
"""
function discover_context_sensitive_call(ast, call_sites :: CallSites, top_level_number, is_top_level, read)
    if typeof(ast) <: Expr
        head = ast.head
        if head == :call || head == :call1
            args = ast.args
            module_name, function_name = resolve_call_names(args)
            arg_types                  = ntuple(i-> type_of_ast_node(args[i+1], call_sites.symbol_info), length(args) - 1)
            item                       = look_for_function(context_sensitive_functions, module_name, function_name, arg_types)
            if item != nothing 
                site = CallSite(ast, Set(), item.fknob_creator, item.fknob_deletor)
                push!(call_sites.sites, site)
                for i in item.matrices
                    push!(site.matrices, typeof(args[i]) == SymbolNode ?
                        args[i].name : args[i])
                end
            end
        end
    end
    return nothing
end

@doc """
Create statements that will create a matrix knob for matrix M.
"""
function create_new_matrix_knob(
    new_stmts :: Vector{Statement},
    M         :: Sym
)
    mknob = gensym(string("mknob", string(M)))
    new_stmt = Expr(:(=), mknob,
                Expr(:call, GlobalRef(SparseAccelerator, :new_matrix_knob)))
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
    fknob = gensym("fknob")
    new_stmt = Expr(:(=), fknob, 
                Expr(:call, GlobalRef(SparseAccelerator, :new_function_knob), 
                    call_site.fknob_creator)
               )
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
    fknob_deletor :: Tuple{Symbol, String},
    fknob         :: Symbol
)
    new_stmt = Expr(:call, GlobalRef(SparseAccelerator, :delete_function_knob), fknob_deletor, fknob)
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
    L         = region.loop
    blocks    = cfg.basic_blocks
    
    # Find call to context-sensitive functions, including their matrix 
    # inputs. Find all definitions of sparse matrices related with the calls.
    call_sites  = CallSites(Set{CallSite}(), symbol_info)
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

            CompilerTools.AstWalker.AstWalk(expr, discover_context_sensitive_call, call_sites)

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
    matrix_knobs         = Dict{Sym, Symbol}()
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
            mknob = matrix_knobs[M]

            # Create statements that will update the knob before every
            # statement that defines the matrix
            for (bb, stmt_idx) in var_defs[M]
                action = InsertBeforeStatement(Vector{Statement}(), bb, stmt_idx)
                push!(actions, action)
                create_increment_matrix_version(action.new_stmts, mknob)
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
    for region in regions
        context_info_discovery(actions, region, func_ast, symbol_info, liveness, cfg)
    end

    dprintln(1, 0, "\nContext-sensitive actions to take:")
    dprintln(1, 1, "", actions)
    
    actions
end

