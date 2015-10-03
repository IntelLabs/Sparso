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
    comment :: AbstractString
end

@doc """"
A pattern not to match with a test's output. The pattern is a regular expression.
The comment explains the pattern, e.g. what the pattern is supposed to do, as 
a friendly message when the pattern does match.
"""
immutable AntiTestPattern
    pattern :: Regex
    comment :: AbstractString
end

immutable Test
    name     :: AbstractString
    command  :: AbstractString
    patterns :: Vector{Union{TestPattern, AntiTestPattern}}
end


root_path       = joinpath(dirname(@__FILE__), "../..")
load_path       = joinpath(root_path, "deps")
julia_command   = joinpath(root_path, "deps", "julia")

# The following patterns often contains "(.|n)*?something". This is to let
# pcregrep (if USE_PCREGREP_REGEX_MATCH is true) match immediately when
# something appears. Otherwise, it will maximize the match and thus (.|n)* can
# match the whole file, causing segmentation fault, etc.

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
        TestPattern(r"New AST:(.|\n)*?return A::[^\s\\]*SparseMatrixCSC.*x::Array",
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

const CG_MATRIX = "bcsstk14.mtx"

const context_test1 = Test(
    "context-test1",
    string("context-test1.jl ", CG_MATRIX),
    [
        TestPattern(r"Original k=208",
                     "Test iterations"
        ),
        TestPattern(r"Original rel_err=8.89\d*e-8",
                     "Test rel_err"
        ),
        TestPattern(r"Constant structures discovered:\n.*\[:A,:L,:U\]",
                     "Test constant structures"
        ),
        TestPattern(r"Structure symmetry discovered:\n.*\[:A\]",
                     "Test structure symmetry"
        ),
        TestPattern(r"L is lower of A",
                     "Test L and A"
        ),
        TestPattern(r"U is upper of A",
                     "Test U and A"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobA.* = \(SparseAccelerator.new_matrix_knob\)\(A,true,true,true,true,false,false\)",
                     "Test if mknobA is generated and is constant valued and constant structured"
        ),
        TestPattern(r"New AST:(.|\n)*?add_mknob_to_fknob\)\(.*mknobA.*,..*fknob.*\)",
                     "Test if mknobA is added to a function knob (for SpMV)"
        ),
        TestPattern(r"New AST:(.|\n)*?Ap = .*:SpMV\)\)\(1,A.*,p.*,.*fknob.*\)",
                     "Test if Ap = A * p has been replaced with SpMV with context info"
        ),
        TestPattern(r"New AST:(.|\n)*?set_reordering_decision_maker",
                     "Test reordering"
        ),
        TestPattern(r"SparseAccelerator.reordering\)\(##fknob#\d*?,##reordering_status#\d*?,(U,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,|A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,){2,2}:__delimitor__,(r,SparseAccelerator.ROW_PERM,?|p,SparseAccelerator.ROW_PERM,?|x,SparseAccelerator.ROW_PERM,?){3,3}\)",
                     "Test reordering"
        ),      
        TestPattern(r"reverse_reordering\)\(##reordering_status#\d*?,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
                     "Test reordering"
        ),
        TestPattern(r"Accelerated k=208",
                     "Test iterations"
        ),
        TestPattern(r"Accelerated rel_err=8.89\d*e-8",
                     "Test rel_err"
        ),
        exception_pattern
    ]
)

const context_test2 = Test(
    "context-test2",
    string("context-test2.jl ", CG_MATRIX, " ", CG_MATRIX),
    [
        TestPattern(r"Original k=208",
                     "Test iterations"
        ),
        TestPattern(r"Original rel_err=8.89\d*e-8",
                     "Test rel_err"
        ),
        TestPattern(r"Manual_context_no_reorder sum of x=40.11\d*",
                     "Test sum of pcg_symgs"
        ),
        TestPattern(r"Manual_context_no_reorder k=208",
                     "Test iterations"
        ),
        TestPattern(r"Manual_context_no_reorder rel_err=8.89\d*e-8",
                     "Test rel_err"
        ),
        TestPattern(r"Manual_context k=208",
                     "Test iterations"
        ),
        TestPattern(r"Manual_context rel_err=8.89\d*e-8",
                     "Test rel_err"
        ),
        exception_pattern
    ]
)

const context_test2_ilu0 = Test(
    "context-test2-ilu0",
    string("context-test2-ilu0.jl  lap3d_4x4x4.mtx"),
    [
        TestPattern(r"Original:(.|\n)*?k=5(.|\n)*?rel_err=4.398690\d*e-9",
                     "Test pcg_symgs"
        ),
        TestPattern(r"With manual context-sensitive optimization without reordering:(.|\n)*?k=5(.|\n)*?rel_err=4.398690\d*e-9",
                     "Test pcg_symgs with contextopt but without reordering"
        ),
        TestPattern(r"With manual context-sensitive optimization:(.|\n)*?k=5(.|\n)*?rel_err=4.398690\d*e-9",
                     "Test pcg_symgs_with_context_opt"
        ),
        exception_pattern
    ]
)

const context_test2_without_reordering = Test(
    "context-test2-without-reordering",
    string("context-test2-without-reordering.jl ", CG_MATRIX),
    [
        TestPattern(r"Original:(.|\n)*?sum of x=40.11\d*(.|\n)*?rel_err=8.89\d*e-8",
                     "Test pcg_symgs"
        ),
        TestPattern(r"With manual context-sensitive optimization:(.|\n)*?sum of x=40.11\d*(.|\n)*?rel_err=8.89\d*e-8",
                     "Test pcg_symgs_with_context_opt"
        ),
        exception_pattern
    ]
)

const context_test3 = Test(
    "context-test3",
    "context-test3.jl ipm/mps/osa-14",
    [
        TestPattern(r"New AST:(.|\n)*?__AT__ = \(Main.ctranspose\)\(A",
                     "Test if accelerated ipm-ref generates AT"
        ),
        TestPattern(r"New AST:(.|\n)*?mknob__AT__.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for AT"
        ),
#        TestPattern(r"New AST:(.|\n)*?mknob__AT__.*new_matrix_knob\)\(__AT__",
#                     "Test if accelerated ipm-ref generates matrix knob for AT"
#        ),
        TestPattern(r"New AST:(.|\n)*?mknobD.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for D"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobA.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for A"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobB.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for B"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobR.*new_matrix_knob",
                     "Test if accelerated ipm-ref generates matrix knob for R"
        ),
        TestPattern(r"New AST:(.|\n)*?new_function_knob(.|\n)*?new_function_knob(.|\n)*?new_function_knob",
                     "Test if accelerated ipm-ref generates function knobs"
        ),
        TestPattern(r"New AST:(.|\n)*?B = .*\(SparseAccelerator,:ADB\)\)\(A.*,D.*,__AT__.*,##fknob#",
                     "Test if accelerated ipm-ref generates ADB"
        ),
        TestPattern(r"New AST:(.|\n)*?B = .*\(SparseAccelerator,:ADB\).*\n.*propagate_matrix_info.*mknobB.*mknobExpr",
                     "Test if accelerated ipm-ref generates propagate_matrix_info after B = ADB"
        ),
        TestPattern(r"New AST:(.|\n)*?R = .*\(SparseAccelerator,:cholfact_int32\)\)\(B.*,##fknob#",
                     "Test if accelerated ipm-ref generates cholfact_int32"
        ),
        TestPattern(r"New AST:(.|\n)*?dy = .*\(SparseAccelerator,:cholfact_inverse_divide\)\)\(R.*,t2.*,##fknob#",
                     "Test if accelerated ipm-ref generates cholfact_inverse_divide"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.set_derivative\)\(.*mknobA.*,SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE,.*mknob__AT.*\)",
                     "Test if accelerated ipm-ref generates a derivative between A and A'"
        ),
        # TODO: enable these patterns
#        TestPattern(r"New AST:(.|\n)*?Rd = .*\(SparseAccelerator,:SpMV\)\)\(1,__AT__,y,1,s-p,fknob_spmv.*\)",
#                     "Test call replacement. Currently cannot match as we need to handle AT specially."
#        ),
        TestPattern(r"New AST:(.|\n)*?Rp = .*?\(SparseAccelerator,:SpMV\)\)\(1,A.*?,x.*?,-1,b.*?,0,.*?fknob.*?\)",
                     "Test call replacement."
        ),
        TestPattern(r"New AST:(.|\n)*?Rc = .*?SparseAccelerator,:element_wise_multiply\)\)\(x.*?,s.*?\)",
                     "Test call replacement."
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPB!\)\)\(Rc,1,Rc,-\(\(Main.min\)\(0.1,100 \* mu",
                     "Test call replacement."
        ),
        TestPattern(r"Original sum of x=715375.98850000",
                     "Test original ipm-ref"
        ),
        TestPattern(r"Original iter 26, resid =  3.93e-15",
                     "Test original ipm-ref"
        ),
        TestPattern(r"Accelerated sum of x=715375.98850000",
                     "Test ipm-ref with context-sensitve optimization"
        ),
        TestPattern(r"Accelerated iter 26, resid =  3.93e-15",
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

const context_test5 = Test(
    "context-test5",
    "context-test5.jl lap3d_4x4x4.mtx",
    [
        TestPattern(r"Original k=5",
                     "Test original pcg_symgs_ilu0 "
        ),
        TestPattern(r"Original rel_err=4.398690\d*e-9",
                     "Test original pcg_symgs_ilu0 "
        ),
        TestPattern(r"Accelerated k=5",
                     "Test iterations"
        ),
        TestPattern(r"Accelerated rel_err=4.398690\d*e-9",
                     "Test rel_err"
        ),
        TestPattern(r"Constant structures discovered:\n.*\[:A,:L,:U\]",
                     "Test constant structures"
        ),
        TestPattern(r"Structure symmetry discovered:\n.*\[:A\]",
                     "Test structure symmetry"
        ),
        TestPattern(r"L is lower of A",
                     "Test L and A"
        ),
        TestPattern(r"U is upper of A",
                     "Test U and A"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobA.* = \(SparseAccelerator.new_matrix_knob\)\(A,true,true,true,true,false,false\)",
                     "Test if mknobA is generated and is constant valued, constant structured, symmetric valued and structured"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobL.* = \(SparseAccelerator.new_matrix_knob\)\(L,true,true,false,false,false,false\)",
                     "Test if mknobL is generated and is constant valued and constant structured"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobU.* = \(SparseAccelerator.new_matrix_knob\)\(U,true,true,false,false,false,false\)",
                     "Test if mknobL is generated and is constant valued and constant structured"
        ),
        TestPattern(r"New AST:(.|\n)*?set_derivative\)\(.*?mknobL.*?,SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC,.*?mknobA.*?\)",
                    "Test if L is found as part of A"
        ),
        TestPattern(r"New AST:(.|\n)*?set_derivative\)\(.*?mknobU.*?,SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC,.*?mknobA.*?\)",
                    "Test if U is found as part of A"
        ),
        TestPattern(r"New AST:(.|\n)*?add_mknob_to_fknob\)\(.*mknobA.*,..*fknob.*\)",
                     "Test if mknobA is added to a function knob (for SpMV)"
        ),
        TestPattern(r"New AST:(.|\n)*?Ap = .*:SpMV\)\)\(1,A.*,p.*,.*fknob.*\)",
                     "Test if Ap = A * p has been replaced with SpMV with context info"
        ),
        TestPattern(r"New AST:(.|\n)*?set_reordering_decision_maker",
                     "Test reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.reordering\)\(##fknob#\d*?,##reordering_status#\d*?,(U,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,|A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,){2,2}:__delimitor__,(r,SparseAccelerator.ROW_PERM,?|p,SparseAccelerator.ROW_PERM,?|x,SparseAccelerator.ROW_PERM,?){3,3}\)",
                     "Test reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?reverse_reordering\)\(##reordering_status#\d*?,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
                     "Test reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?r = b.* - .*SparseAccelerator,:SpMV\)\)\(A.*,x.*\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?normr0 = .*SparseAccelerator,:norm\)\)\(r.*\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?rz = .*SparseAccelerator,:dot\)\)\(r.*,z.*\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?alpha = old_rz.*?/.*?SparseAccelerator,:dot\)\)\(p.*?,Ap.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY!\)\)\(x,1,x.*?,alpha.*?,p.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY!\)\)\(r,1,r.*?,-alpha.*?,Ap.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?rel_err = .*?SparseAccelerator,:norm\)\)\(r.*?\).*?/ normr0",
                     "Test call replacement"
        ),
# TODO: enable this test.
#        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:copy!\)\)\(z.*?,r.*?\)",
#                     "Test call replacement"
#        )
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:fwdTriSolve!\)\)\(L.*?,z.*?,.*?fknob.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:bwdTriSolve!\)\)\(U.*?,z.*?,.*?fknob.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY!\)\)\(p,1,z.*?,beta.*?,p.*?\)",
                     "Test call replacement"
        ),
        exception_pattern
    ]
)

const pagerank_test1 = Test(
    "pagerank-test1",
    "pagerank.jl  hmatrix.1024.mtx",
    [
        TestPattern(r"Original:(.\n)*?\s*error = 1.54\d*e-8",
                     "Test original pagerank"
        ),
        TestPattern(r"Accelerated:(.\n)*?\s*error = 1.54\d*e-8",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?set_reordering_decision_maker",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?Ap = \(\(top\(getfield\)\)\(SparseAccelerator,:SpMV\)\)\(1 - r.*?,A.*?,p.*?,0,p.*?,r.*?,##fknob.*?\)",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"reverse_reordering\)\(##reordering_status#\d*?,:__delimitor__,p,SparseAccelerator.ROW_PERM\)",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"err = .*?\(SparseAccelerator,:norm\)\)\(Ap.*? - p.*?\).*? / .*?\(SparseAccelerator,:norm\)\)\(p.*?\)",
                     "Test pagerank if call replacement of norm happens"
        ),
        exception_pattern
    ]
)

const cosp2_test1 = Test(
    "cosp2-test1",
    "cosp2.jl",
    [
        TestPattern(r"Original:\nD Sparsity AAN = 27050.134369218962, fraction = 0.0001791459611337646 avg = 2.2013455704116995, max = 2.0122570097076777\nNumber of iterations = 39(.|\n)*?End original.",
                     "Test original"
        ),
        TestPattern(r"CoSP2_call_replacement:\nD Sparsity AAN = 12212.785128790038, fraction = 8.088207992441149e-5 avg = 0.9938789981111684, max = 1.2808837088549982\nNumber of iterations = 25(.|\n)*?End CoSP2_call_replacement.",
                     "Test CoSP2_call_replacement"
        ),
        TestPattern(r"CoSP2_call_replacement_and_context_opt:\nD Sparsity AAN = 12212.785128790038, fraction = 8.088207992441149e-5 avg = 0.9938789981111684, max = 1.2808837088549991\nNumber of iterations = 25(.|\n)*?End CoSP2_call_replacement_and_context_opt",
                     "Test CoSP2_call_replacement_and_context_opt"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobX.*? = \(SparseAccelerator.new_matrix_knob\)\(false,.*?,true,true,false,false\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.add_mknob_to_fknob\)\(.*?mknobX.*?,.*?fknob.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?trX = \(\(top\(getfield\)\)\(SparseAccelerator,:trace\)\)\(X.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?X2 = \(\(top\(getfield\)\)\(SparseAccelerator,:SpSquareWithEps\)\)\(X.*?,eps.*?,.*?fknob.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?trX2 = \(\(top\(getfield\)\)\(SparseAccelerator,:trace\)\)\(X2.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?X = \(\(top\(getfield\)\)\(SparseAccelerator,:SpAdd\)\)\(2,X.*?,-1,X2.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?D Sparsity AAN = 12212.785128790038, fraction = 8.088207992441149e-5 avg = 0.9938789981111684, max = 1.2808837088549991\nNumber of iterations = 25(.|\n)*?End accelerated",
                     "Test accelerated"
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
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.*bigM.* A .*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets bigM"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.* A .*bigM.*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets bigM"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.*Rd.* A .*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets Rd"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.* A .*Rd.*\) Use\(",
                        "Test liveness for ipm-ref: A should not be updated in the block that sets Rd"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.* mu .*\) Use\(.*\n.* = mu <=",
                        "Test liveness for ipm-ref: mu should not be updated in the block that tests mu <= 1.0e-7"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.*blas1_time.* relResidual .*\) Use\(.*\n\s*blas1_time =.*\n.*Ac_mul_B",
                        "Test liveness for ipm-ref: relResidual, x, p should not be updated in the block that sets blas1_time"
        ),
        AntiTestPattern(r"Liveness of basic blocks:(.|\n)*?Def\(.* relResidual .*blas1_time.*\) Use\(.*\n\s*blas1_time =.*\n.*Ac_mul_B",
                        "Test liveness for ipm-ref: relResidual, x, p should not be updated in the block that sets blas1_time"
        ),
        exception_pattern
    ]
)

const call_replacement_test1 = Test(
    "call-replacement-test1",
    "call-replacement-test1.jl small-diag.mtx",
    [
        TestPattern(r"AST:(.|\n)*?Main.dot(.|\n)*?New AST(.|\n)*?SparseAccelerator,:dot",
                     "Test call replacement of Main.dot with SparseAccelerator.dot."
        ),
        exception_pattern
    ]
)

const call_replacement_test2 = Test(
    "call-replacement-test2",
    "call-replacement-test2.jl small-diag.mtx",
    [
        TestPattern(r"AST:(.|\n)*?A::[^\s\\]*SparseMatrixCSC.*\* x::Array(.|\n)*?New AST(.|\n)*?SparseAccelerator,:SpMV.*A::[^\s\\]*SparseMatrixCSC.*,x::Array",
                     "Test call replacement of * with SparseAccelerator.SpMV."
        ),
        exception_pattern
    ]
)

const call_replacement_test3 = Test(
    "call-replacement-test3",
    "call-replacement-test3.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:SpMV\!\)\)\(y::Array\{Float64,1\},A::[^\s\\]*SparseMatrixCSC\{Float64,Int64\},x::Array\{Float64,1\}\)",
                     "Test call replacement of SpMV! for A_mul_B!(y, A, x)."
        ),
        exception_pattern
    ]
)

const call_replacement_test4 = Test(
    "call-replacement-test4",
    "call-replacement-test4.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:SpMV\!\)\)\(y::Array\{Float64,1\},0.1,A::[^\s\\]*SparseMatrixCSC\{Float64,Int64\},x::Array\{Float64,1\},0.1,y::Array\{Float64,1\},0.0\)",
                     "Test call replacement of SpMV for A_mul_B!(0.1, A, x, 0.1, y)."
        ),
        exception_pattern
    ]
)

const call_replacement_test5 = Test(
    "call-replacement-test5",
    "call-replacement-test5.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY\!\)\)\(x,1,x::Array\{Float64,1\},alpha::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for x += alpha * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test6 = Test(
    "call-replacement-test6",
    "call-replacement-test6.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY\!\)\)\(x,1,x::Array\{Float64,1\},-alpha::Float64::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for x -= alpha * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test7 = Test(
    "call-replacement-test7",
    "call-replacement-test7.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY\!\)\)\(p,1,r::Array\{Float64,1\},beta::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for p = r + beta * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test8 = Test(
    "call-replacement-test8",
    "call-replacement-test8.jl small-diag.mtx",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:SpMV\!\)\)\(p,1 - r::Float64::Float64,A::[^\s\\]*SparseMatrixCSC\{Float64,Int32\},p::Array\{Float64,1\},0,p::Array\{Float64,1\},r::Float64\)",
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
        TestPattern(r"sum of x=-1.5773120434107\d*e-5",
                     "Test orig sum"
        ),

        TestPattern(r"accel sum of x=-1.5773120434107\d*e-5",
                     "Test accelerated sum"
        ),
        exception_pattern
    ]
)

const call_replacement_test11 = Test(
    "call-replacement-test11",
    "call-replacement-test11.jl small-diag.mtx",
    [
        TestPattern(r"sum of x=-1.5773120434107\d*e-5",
                     "Test orig sum"
        ),

        TestPattern(r"accel sum of x=-1.5773120434107\d*e-5",
                     "Test accelerated sum"
        ),
        exception_pattern
    ]
)

const call_replacement_test12 = Test(
    "call-replacement-test12",
    string("call-replacement-test12.jl  ", CG_MATRIX),
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:WAXPBY\!\)\)\(p,1,z::Array\{Float64,1\},beta::Float64,p::Array\{Float64,1\}\)",
                     "Test call replacement of WAXPBY! for p = r + beta * p."
        ),
        exception_pattern
    ]
)

const name_resolution_test1 = Test(
    "name-resolution-test1",
    "name-resolution-test1.jl small-diag.mtx",
    [
        TestPattern(r"Module name: X\.Y\.Z\.U\.V\.W\nFunction name: f(.|\n)*?Module name: Main\nFunction name: \*",
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
    "constant-value-test1.jl ipm/mps/osa-14",
    [
        TestPattern(r"Constants discovered:.*\n.*\[(:A,?|:b,?|:p,?){3,3}\]",
                     "Test ipm-ref that A, b and p are recognized as loop constants."
        ),
        exception_pattern
    ]
)

const single_def_test1 = Test(
    "single-def-test1",
    "single-def-test1.jl",
    [
        TestPattern(Regex("Single-defs discovered:.*\\n.*\[.*[^:]*:B[^:].*\]"),
                     "Test ipm-ref that B is recognized as single-defs in the loop."
        ),
        TestPattern(Regex("Single-defs discovered:.*\\n.*\[.*[^:]*:D[^:].*\]"),
                     "Test ipm-ref that D is recognized as single-defs in the loop."
        ),
        TestPattern(Regex("Single-defs discovered:.*\\n.*\[.*[^:]*:R[^:].*\]"),
                     "Test ipm-ref that R is recognized as single-defs in the loop."
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

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
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
        TestPattern(Regex("Func Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Func Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Func Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that R is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
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
    "constant-structure-test1.jl ipm/mps/osa-14",
    [
        TestPattern(Regex("Loop-4 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B, :D, :R])),
                     "Test ipm-ref that A B D R are recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const symmetric_value_test1 = Test(
    "symmetric-value-test1",
    "symmetric-value-test1.jl",
    [
        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:B, :E, :F, :G, :S])),
                     "Test ipm-ref that B E F G S is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const symmetric_value_test2 = Test(
    "symmetric-value-test2",
    "symmetric-value-test2.jl",
    [
        TestPattern(r"Value symmetry discovered:.*\n.*\[.*:A.*\]",
                     "Test ipm-ref that A is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const symmetric_value_test3 = Test(
    "symmetric-value-test3",
    "symmetric-value-test3.jl",
    [
        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:B, :C, :G, :S])),
                     "Test ipm-ref that B E F G S is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const symmetric_value_test4 = Test(
    "symmetric-value-test4",
    "symmetric-value-test4.jl",
    [
        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test ipm-ref that B E F G S is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)


const symmetric_structure_test1 = Test(
    "symmetric-structure-test1",
    "symmetric-structure-test1.jl",
    [
        TestPattern(Regex("Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as symmetric in structure."
        ),
        exception_pattern
    ]
)

const lower_upper_test1 = Test(
    "lower-upper-test1",
    "lower-upper-test1.jl",
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
    context_test2_ilu0,
    context_test3,
    context_test4,
    context_test5,
    pagerank_test1,
    cosp2_test1,
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
    single_def_test1,
    constant_value_test1,
    set_matrix_property_test1,
    set_matrix_property_test2,
    set_matrix_property_test3,
    set_matrix_property_test4,
    set_matrix_property_test5,
    constant_structure_test1,
    symmetric_value_test1,
    symmetric_value_test2,
    symmetric_value_test3,
    symmetric_value_test4,
#    symmetric_structure_test1,
    lower_upper_test1,
]

const fast_tests = [
#    context_test1,
    context_test2,
    context_test2_without_reordering,
    context_test3,
    context_test4,
    context_test5,
    pagerank_test1,
    cosp2_test1
]

# If true, use pcregrep for regular expression match. 
# If false, use Julia (If PCRE JIT stack overflow: enlarge JIT_STACK_MAX_SIZE in
# julia/base/pcre.jl (times it with 10) and retry).
const USE_PCREGREP_REGEX_MATCH = true

# Run tests with multiple threads?
const USE_THREADS = true

function get_julia_ver()
    s, p = open(`$julia_command -v`)
    readline(s)
end

if !isreadable(julia_command)
    error("Please install (softlink) julia command to \"" * julia_command  * "\".")
elseif !ismatch(r"\.*0.4.0-rc3", get_julia_ver())
    error("Wrong julia version! 0.4.0-rc3 is required!")
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
    else
        matched_tests = []
        for t in all_tests
            if startswith(t.name, ARGS[1])
                push!(matched_tests, t)
            end
        end
        if !isempty(matched_tests)
            tests = matched_tests
        else
            print("No testing group matched.\n")
        end
    end  
    print(length(tests), " test(s).\n")
end

ENV["JULIA_LOAD_PATH"] = root_path * "/deps"
# avoid package precompilation issue
run(`$julia_command precompl.jl`)

@doc """
Run a test
"""
function run_test(
    test :: Test
)
    log  = string(test.name, ".log")
    
    # Run the command. Redirect output to the log file
    split_res = split(test.command)
    cmd = `$julia_command $split_res`
    successful = success(pipeline(cmd, stdout=log, stderr=log, append=false))
   
    # Some patterns are not intended to run, but to check compiler outputs
    # expected code
    # TODO: maybe we should add a field to a test to indicate if it should be
    # executed?
    successful = true

    if successful
        # Match with patterns
        log_file = open(log, "a")
        for pattern in test.patterns
            assert(typeof(pattern) == TestPattern || typeof(pattern) == AntiTestPattern)
            if !USE_PCREGREP_REGEX_MATCH
                # Read the output to a string
                output = open(readall, log)
                m       = match(pattern.pattern, output)
                matched = (m != nothing)
            else
                pattern_str = pattern.pattern.pattern

                cmd  = `pcregrep --buffer-size=10000000 -M $pattern_str $log`
                matched = success(pipeline(cmd, stdout=DevNull, stderr=DevNull))
            end
            if (!matched && typeof(pattern) == TestPattern) ||
               ( matched && typeof(pattern) == AntiTestPattern)
                comment = pattern.comment
                write(log_file, "\n****** Failed in ", 
                    (typeof(pattern) == AntiTestPattern) ? "anti-pattern\n\t" : "pattern\n\t",
                    string(pattern.pattern), "\n\tComment: ", comment)
                successful = false
            end
        end
        close(log_file)
    end

    if successful
        #rm(log)
    end

    successful
end

total = length(tests)
succ = 0

if USE_THREADS == false
    for test in tests
        print("Testing ", test.name)
        s = run_test(test)
        if s
            succ = succ + 1
            println(": Pass")
        else
            println(": FAIL. See ", test.name * ".log")
        end
    end
else
    tasks = []
    task_map = Dict()
    for test in tests
        push!(tasks, @schedule(run_test(test)))
        task_map[last(tasks)] = test
        sleep(0.5)
    end
    finished = []
    running  = []
    passed   =  Dict()
    while length(finished) < total
        for t in tasks
            if istaskdone(t)
                if !in(t, finished)
                    push!(finished, t)
                    name = task_map[t].name
                    ret = wait(t)
                    print(length(finished), "/", total, " ")
                    if ret 
                        println(name, ": Pass")
                        passed[name] = true
                    else
                        println(name, ": FAIL. See ", name * ".log")
                        passed[name] = false
                    end
                    succ = succ + ret
                end
            elseif istaskstarted(t)
                if !in(t, running)
                    push!(running, t)
                    name = task_map[t].name
                    println("Testing ", name)
                end
            end
        end
        sleep(0.5)
    end

    # Print out the results in a fixed order for easier checking
    println("\nFinal report:")
    for test in tests
        println(test.name, ": ", passed[test.name] ? "Pass" : "FAIL. See $(test.name).log")
    end
end

println("Total: ", total)
println("Pass : ", succ)
println("Fail : ", total-succ)

