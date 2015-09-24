# This file is used for preventing regressions and keeping the quality of the 
# repository higher every checkin. 
# For the module you write, please write as many unit test cases as possible, and
# add the corresponding entries for them -- so that your code will not be accidently
# affected by other people's checkins!

@doc """"
A pattern to match with a test's output. The pattern is a regular expression.
The comment explains the pattern, e.g. what the pattern is supposed to do, as 
a friendly message when the pattern does not match.
"""
immutable TestPattern
    pattern :: Regex
    comment :: String
end

@doc """"
A pattern not to match with a test's output. The pattern is a regular expression.
The comment explains the pattern, e.g. what the pattern is supposed to do, as 
a friendly message when the pattern does match.
"""
immutable AntiTestPattern
    pattern :: Regex
    comment :: String
end

immutable Test
    name     :: String
    command  :: String
    patterns :: Vector{Union(TestPattern, AntiTestPattern)}
end


root_path       = joinpath(dirname(@__FILE__), "../..")
load_path       = joinpath(root_path, "dpes")
julia_command   = joinpath(root_path, "deps", "julia")

const exception_pattern = AntiTestPattern(
    r"Exception!",
    "Test if an exception has been thrown in sparse accelerator"
)

const sanity_test1 = Test(
    "sanity-test1",
    "sanity-test1.jl",
    [
        TestPattern(r" {4,4}1 tab",
                     "Test tabbed output facility: 1 tab."
        ),
        TestPattern(r" {8,8}2 tabs",
                     "Test tabbed output facility: 2 tabs."
        ),
        TestPattern(r" {12,12}3 tabs",
                     "Test tabbed output facility: 3 tabs."
        ),
        exception_pattern
    ]
)

const sanity_test2 = Test(
    "sanity-test2",
    "sanity-test2.jl",
    [
        TestPattern(r"I do nothing!",
                     "Test if optimization framework adds a new optimization pass."
        ),
        exception_pattern
    ]
)

const sanity_test3 = Test(
    "sanity-test3",
    "sanity-test3.jl",
    [
        TestPattern(r"New AST:(.|\n)*return A::[^\s\\]*SparseMatrixCSC.*x::Array",
                     "Test if optimization framework invokes SparseAccelerator to generate a new AST."
        ),
        exception_pattern
    ]
)

const spmv_sanity_test1 = Test(
    "spmv-sanity-test1",
    "spmv-sanity-test1.jl tiny-diag.mtx",
    [
        TestPattern(r"Original sum of p=20.169999999999995",
                     "Test original code"
        ),
        TestPattern(r"Manual sum of p=20.169999999999995",
                     "Test code manually replaced with SpMV!"
        ),
        exception_pattern
    ]
)

const context_test1 = Test(
    "context-test1",
    "context-test1.jl small-diag.mtx",
    [
        TestPattern(r"Original sum of x=-1.57731204341073\d*e-5",
                     "Test sum of pcg_symgs"
        ),
        TestPattern(r"Original k=4",
                     "Test iterations"
        ),
        TestPattern(r"Original rel_err=6.38\d*e-13",
                     "Test rel_err"
        ),
        TestPattern(r"New AST:(.|\n)*mknobA.* = \(SparseAccelerator.new_matrix_knob\)\(A,true,true,false,false,false,false\)",
                     "Test if mknobA is generated and is constant valued and constant structured"
        ),
        TestPattern(r"New AST:(.|\n)*add_mknob_to_fknob\)\(.*mknobA.*,..*fknob.*\)",
                     "Test if mknobA is added to a function knob (for SpMV)"
        ),
        TestPattern(r"New AST:(.|\n)*Ap = .*:SpMV\)\)\(A.*,p.*,.*fknob.*\)",
                     "Test if Ap = A * p has been replaced with SpMV with context info"
        ),
        TestPattern(r"Accelerated sum of x=-1.577312043\d*e-5",
                     "Test sum of pcg_symgs"
        ),
        TestPattern(r"Accelerated k=4",
                     "Test iterations"
        ),
        TestPattern(r"Accelerated rel_err=6.38\d*e-13",
                     "Test rel_err"
        ),
        exception_pattern
    ]
)

const context_test2 = Test(
    "context-test2",
    "context-test2.jl small-diag.mtx small-diag.mtx",
    [
        TestPattern(r"Original sum of x=-1.577312043\d*e-5",
                     "Test sum of pcg_symgs"
        ),
        TestPattern(r"Original k=4",
                     "Test iterations"
        ),
        TestPattern(r"Original rel_err=6.38\d*e-13",
                     "Test rel_err"
        ),
        TestPattern(r"Manual_context_no_reorder sum of x=-1.577312043\d*e-5",
                     "Test sum of pcg_symgs"
        ),
        TestPattern(r"Manual_context_no_reorder k=4",
                     "Test iterations"
        ),
        TestPattern(r"Manual_context_no_reorder rel_err=\d*e-11",
                     "Test rel_err"
        ),
        TestPattern(r"Manual_context sum of x=-1.577312043\d*e-5",
                     "Test sum of pcg_symgs"
        ),
        TestPattern(r"Manual_context k=4",
                     "Test iterations"
        ),
        TestPattern(r"Manual_context rel_err=\d*e-11",
                     "Test rel_err"
        ),
        exception_pattern
    ]
)

const context_test2_without_reordering = Test(
    "context-test2-without-reordering",
    "context-test2-without-reordering.jl small-diag.mtx",
    [
        TestPattern(r"Original:(.|\n)*sum of x=-1.5773120434\d*e-5(.|\n)*rel_err=\d.\d*e-13",
                     "Test pcg_symgs"
        ),
        TestPattern(r"With manual context-sensitive optimization:(.|\n)*sum of x=-1.5773120434\d*e-5(.|\n)*rel_err=\d.\d*e-11",
                     "Test pcg_symgs_with_context_opt"
        ),
        exception_pattern
    ]
)

const context_test3 = Test(
    "context-test3",
    "context-test3.jl ipm/mps/osa-14",
    [
        TestPattern(r"New AST:(.|\n)*__AT__ = \(Main.ctranspose\)\(A",
                     "Test if accelerated ipm-ref generates AT"
        ),
        TestPattern(r"New AST:(.|\n)*mknob__AT__.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for AT"
        ),
#        TestPattern(r"New AST:(.|\n)*mknob__AT__.*new_matrix_knob\)\(__AT__",
#                     "Test if accelerated ipm-ref generates matrix knob for AT"
#        ),
        TestPattern(r"New AST:(.|\n)*mknobD.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for D"
        ),
        TestPattern(r"New AST:(.|\n)*mknobA.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for A"
        ),
        TestPattern(r"New AST:(.|\n)*mknobB.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for B"
        ),
        TestPattern(r"New AST:(.|\n)*mknobR.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for R"
        ),
        TestPattern(r"New AST:(.|\n)*new_function_knob(.|\n)*new_function_knob(.|\n)*new_function_knob",
                     "Test if accelerated ipm-ref generates function knobs"
        ),
        TestPattern(r"New AST:(.|\n)*B = .*\(SparseAccelerator,:ADB\)\)\(A.*,D.*,__AT__.*,##fknob#",
                     "Test if accelerated ipm-ref generates ADB"
        ),
        TestPattern(r"New AST:(.|\n)*B = .*\(SparseAccelerator,:ADB\).*\n.*propagate_matrix_info.*mknobB.*mknobExpr",
                     "Test if accelerated ipm-ref generates propagate_matrix_info after B = ADB"
        ),
        TestPattern(r"New AST:(.|\n)*R = .*\(SparseAccelerator,:cholfact_int32\)\)\(B.*,##fknob#",
                     "Test if accelerated ipm-ref generates cholfact_int32"
        ),
        TestPattern(r"New AST:(.|\n)*dy = .*\(SparseAccelerator,:cholfact_inverse_divide\)\)\(R.*,t2,##fknob#",
                     "Test if accelerated ipm-ref generates cholfact_inverse_divide"
        ),
        TestPattern(r"Original sum of x=715375.9885000014",
                     "Test original ipm-ref"
        ),
        TestPattern(r"Accelerated sum of x=715375.9885000014",
                     "Test ipm-ref with context-sensitve optimization"
        ),
        exception_pattern
    ]
)

const context_test4 = Test(
    "context-test4",
    "context-test4.jl ipm/mps/osa-14",
    [
        TestPattern(r"Original sum of x=715375.988500001",
                     "Test original ipm-ref"
        ),
        TestPattern(r"Manual_context sum of x=715375.988500001",
                     "Test ipm-ref with context-sensitve optimization"
        ),
        exception_pattern
    ]
)

const liveness_test1 = Test(
    "liveness-test1",
    "liveness-test1.jl small-diag.mtx",
    [
        TestPattern(r"Def",
                     "Test liveness for cg."
        ),
        exception_pattern
    ]
)

const liveness_test2 = Test(
    "liveness-test2",
    "liveness-test2.jl",
    [
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.*bigM.* A .*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets bigM"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.* A .*bigM.*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets bigM"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.*Rd.* A .*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets Rd"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.* A .*Rd.*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets Rd"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.* mu .*\) Use\(.*\n.* = mu <=",
                        "Test liveness for ipm-ref: mu should not be updated in the block that tests mu <= 1.0e-7"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.*blas1_time.* relResidual .*\) Use\(.*\n\s*blas1_time =.*\n.*Ac_mul_B",
                        "Test liveness for ipm-ref: relResidual, x, p should not be updated in the block that sets blas1_time"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*Def\(.* relResidual .*blas1_time.*\) Use\(.*\n\s*blas1_time =.*\n.*Ac_mul_B",
                        "Test liveness for ipm-ref: relResidual, x, p should not be updated in the block that sets blas1_time"
        ),
        exception_pattern
    ]
)

const call_replacement_test1 = Test(
    "call-replacement-test1",
    "call-replacement-test1.jl small-diag.mtx",
    [
        TestPattern(r"AST:(.|\n)*Main.dot(.|\n)*New AST(.|\n)*SparseAccelerator,:dot",
                     "Test call replacement of Main.dot with SparseAccelerator.dot."
        ),
        exception_pattern
    ]
)

const call_replacement_test2 = Test(
    "call-replacement-test2",
    "call-replacement-test2.jl small-diag.mtx",
    [
        TestPattern(r"AST:(.|\n)*A::[^\s\\]*SparseMatrixCSC.*\* x::Array(.|\n)*New AST(.|\n)*SparseAccelerator,:SpMV.*A::[^\s\\]*SparseMatrixCSC.*,x::Array",
                     "Test call replacement of * with SparseAccelerator.SpMV."
        ),
        exception_pattern
    ]
)

const call_replacement_test3 = Test(
    "call-replacement-test3",
    "call-replacement-test3.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:SpMV\!\)\)\(y::Array\{Float64,1\},A::[^\s\\]*SparseMatrixCSC\{Float64,Int64\},x::Array\{Float64,1\}\)",
                     "Test call replacement of SpMV! for A_mul_B!(y, A, x)."
        ),
        exception_pattern
    ]
)

const call_replacement_test4 = Test(
    "call-replacement-test4",
    "call-replacement-test4.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:SpMV\!\)\)\(y::Array\{Float64,1\},0.1,A::[^\s\\]*SparseMatrixCSC\{Float64,Int64\},x::Array\{Float64,1\},0.1,y::Array\{Float64,1\},0.0\)",
                     "Test call replacement of SpMV for A_mul_B!(0.1, A, x, 0.1, y)."
        ),
        exception_pattern
    ]
)

const call_replacement_test5 = Test(
    "call-replacement-test5",
    "call-replacement-test5.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:WAXPBY\!\)\)\(x,1,x::Array\{Float64,1\},alpha::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for x += alpha * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test6 = Test(
    "call-replacement-test6",
    "call-replacement-test6.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:WAXPBY\!\)\)\(x,1,x::Array\{Float64,1\},-alpha::Float64::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for x -= alpha * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test7 = Test(
    "call-replacement-test7",
    "call-replacement-test7.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:WAXPBY\!\)\)\(p,1,r::Array\{Float64,1\},beta::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for p = r + beta * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test8 = Test(
    "call-replacement-test8",
    "call-replacement-test8.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:SpMV\!\)\)\(p,1 - r::Float64::Float64,A::[^\s\\]*SparseMatrixCSC\{Float64,Int32\},p::Array\{Float64,1\},0,p::Array\{Float64,1\},r::Float64\)",
                     "Test call replacement of SpMV! in simple page rank."
        ),
        exception_pattern
    ]
)

const call_replacement_test9 = Test(
    "call-replacement-test9",
    "call-replacement-test9.jl small-diag.mtx",
    [
        TestPattern(r"sum of x=-1.57731204341073\d*e-5",
                     "Test orig sum"
        ),

        TestPattern(r"accel sum of x=-1.57731204341073\d*e-5",
                     "Test accelerated sum"
        ),
        
        exception_pattern
    ]
)

const call_replacement_test10 = Test(
    "call-replacement-test10",
    "call-replacement-test10.jl small-diag.mtx",
    [
        TestPattern(r"sum of x=-1.5773120434107317e-5",
                     "Test orig sum"
        ),

        TestPattern(r"accel sum of x=-1.57731204341073\d*e-5",
                     "Test accelerated sum"
        ),
        exception_pattern
    ]
)

const call_replacement_test11 = Test(
    "call-replacement-test11",
    "call-replacement-test11.jl small-diag.mtx",
    [
        TestPattern(r"sum of x=-1.57731204341073\d*e-5",
                     "Test orig sum"
        ),

        TestPattern(r"accel sum of x=-1.57731204341073\d*e-5",
                     "Test accelerated sum"
        ),
        exception_pattern
    ]
)

const call_replacement_test12 = Test(
    "call-replacement-test12",
    "call-replacement-test12.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*SparseAccelerator,:WAXPBY\!\)\)\(p,1,z::Array\{Float64,1\},beta::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for p = r + beta * p."
        ),
        exception_pattern
    ]
)

const name_resolution_test1 = Test(
    "name-resolution-test1",
    "name-resolution-test1.jl small-diag.mtx",
    [
        TestPattern(r"Module name: X\.Y\.Z\.U\.V\.W\nFunction name: f(.|\n)*Module name: Main\nFunction name: \*",
                     "Test name resolution."
        ),
        exception_pattern
    ]
)

#
function gen_set_regex_string(
    set :: Array
) 
    ret = "\\["
    for e in set
        ret = ret * "[^:]*:" * string(e)
    end
    ret = ret * "[^:]*\\]"
    ret 
end


const constant_value_test1 = Test(
    "constant-value-test1",
    "constant-value-test1.jl",
    [
        TestPattern(r"Constants discovered:.*\n.*\[.*:A.*\]",
                     "Test ipm-ref that A is recognized as a loop constant."
        ),
        exception_pattern
    ]
)

const single_def_test1 = Test(
    "single-def-test1",
    "single-def-test1.jl",
    [
        TestPattern(Regex("Single-defs discovered:.*\\n.*" * gen_set_regex_string([:B, :D, :R])),
                     "Test ipm-ref that B D R are recognized as single-defs in the loop."
        ),
        exception_pattern
    ]
)

const set_matrix_property_test1 = Test(
    "set-matrix-property-test1",
    "set-matrix-property-test1.jl",
    [
        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that R is recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const set_matrix_property_test2 = Test(
    "set-matrix-property-test2",
    "set-matrix-property-test2.jl",
    [

        TestPattern(Regex("Loop0 Matrix structures discovered:.*\\n.*" * gen_set_regex_string([:A, :L, :U])),
                     "Test ipm-ref that A L U are recognized as matrics in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :L, :U])),
                     "Test ipm-ref that A L U are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that R is recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const set_matrix_property_test3 = Test(
    "set-matrix-property-test3",
    "set-matrix-property-test3.jl",
    [
        TestPattern(Regex("Func Constant structures discovered:.*\\n.*" * gen_set_regex_string([:B])),
                     "Test ipm-ref that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Func Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Func Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that R is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that R is recognized as constant in structure."
        ),
        exception_pattern
    ]
)


const set_matrix_property_test4 = Test(
    "set-matrix-property-test4",
    "set-matrix-property-test4.jl",
    [

        TestPattern(Regex("Func Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that A is recognized as matrics in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that A is recognized as matrics in structure."
        ),

        TestPattern(Regex("Loop3 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as constant in structure."
        ),

        exception_pattern
    ]
)

const set_matrix_property_test5 = Test(
    "set-matrix-property-test5",
    "set-matrix-property-test5.jl",
    [
        TestPattern(r"Func Upper/Lower matrix discovered:.*\n[\s]*A is lower of B",
                     "Test ipm-ref that A is recognized as matrics in structure."
        ),
        TestPattern(r"Loop0 Upper/Lower matrix discovered:.*\n[\s]*A is lower of B",
                     "Test ipm-ref that A is recognized as matrics in structure."
        ),
        TestPattern(r"Loop3 Upper/Lower matrix discovered:.*\n[\s]*A is lower of B\n[\s]*B is upper of C",
                     "Test ipm-ref that A is recognized as matrics in structure."
        ),
        exception_pattern
    ]
)


const constant_structure_test1 = Test(
    "constant-structure-test1",
    "constant-structure-test1.jl",
    [
        TestPattern(Regex("Loop-4 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B, :D, :R])),
                     "Test ipm-ref that A B D R are recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const value_symmetry_test1 = Test(
    "value-symmetry-test1",
    "value-symmetry-test1.jl",
    [
        TestPattern(r"Value symmetry discovered:.*\n.*\[.*:A.*\]",
                     "Test ipm-ref that A is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const structure_symmetry_test1 = Test(
    "structure-symmetry-test1",
    "structure-symmetry-test1.jl",
    [
        TestPattern(Regex("Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as symmetric in structure."
        ),
        exception_pattern
    ]
)

const all_tests = [
    sanity_test1,
    sanity_test2,
    sanity_test3,
    spmv_sanity_test1,
    context_test1,
    context_test2,
    context_test2_without_reordering,
    context_test3,
    context_test4,
    liveness_test1,
    liveness_test2,
    call_replacement_test1,
    call_replacement_test2,
    call_replacement_test3,
    call_replacement_test4,
    call_replacement_test5,
    call_replacement_test6,
    call_replacement_test7,
    call_replacement_test8,
    call_replacement_test9,
    call_replacement_test10,
    call_replacement_test11,
    call_replacement_test12,
    name_resolution_test1,
    constant_value_test1,
    single_def_test1,
    set_matrix_property_test1,
    set_matrix_property_test2,
    set_matrix_property_test3,
    set_matrix_property_test4,
    set_matrix_property_test5,
    constant_structure_test1,
    value_symmetry_test1,
    structure_symmetry_test1
]

const fast_tests = [
    context_test2,
    context_test2_without_reordering,
    context_test3,
    context_test4,
]

function get_julia_ver()
    s, p = open(`$julia_command -v`)
    readline(s)
end

if !ismatch(r"\.*0.4.0-rc1", get_julia_ver())
    error("Wrong julia version! 0.4.0-rc1 is required!")
end

if isreadable("regression.conf")
    include("regression.conf")
    tests = local_tests
else
    tests = fast_tests
end

if length(ARGS) > 0
    if ARGS[1] == "all"
        tests = all_tests
    elseif ARGS[1] == "fast"
        tests = fast_tests
    elseif ARGS[1] == "none"
        tests = []
    end 
end


ENV["JULIA_LOAD_PATH"] = root_path * "/deps"
# avoid package precompilation issue
run(`$julia_command precompl.jl`)

fail = 0
succ = 0
old_stderr = STDERR
old_stdout = STDOUT
for test in tests
    print("Testing ", test.name)
    log  = string(test.name, ".log")
    
    # Run the command. Redirect output to the log file
    #file = open(log, "w+")
    #redirect_stderr(file)
    #redirect_stdout(file)
    output = ""
    try
        split_res = split(test.command)
        cmd = `$julia_command $split_res`
        #run(cmd |> log)
        run(pipeline(cmd, stdout=log, stderr=log))
    catch ex
        println("\nexception = ", ex)
    finally
        #flush(file)
    end
    #close(file)
    #redirect_stderr(old_stderr)
    #redirect_stdout(old_stdout)

    # Read the output to a string
    output = open(readall, log)
    
    # Match with patterns
    successful = true
    for pattern in test.patterns
        assert(typeof(pattern) == TestPattern || typeof(pattern) == AntiTestPattern)

        # This may cause PCRE JIT stack overflow. So change to use grep
        # m = match(pattern.pattern, output)
        run(`grep -E $pattern.pattern $log`)

        if (m == nothing && typeof(pattern) == TestPattern) ||
           (m != nothing && typeof(pattern) == AntiTestPattern)
            comment = pattern.comment
            file = open(log, "a")
            write(file, "\n****** Failed in ", 
                (typeof(pattern) == AntiTestPattern) ? "anti-pattern\n\t" : "pattern\n\t",
                string(pattern.pattern), "\n\tComment: ", comment)
            close(file)
            successful = false
        end
    end
    if successful
        succ = succ + 1
        rm(log)
        println(": Pass")
    else
        fail = fail + 1
        println(": FAIL. See ", log)
    end
end

println("Total: ", fail + succ)
println("Pass : ", succ)
println("Fail : ", fail)
flush(STDOUT)
