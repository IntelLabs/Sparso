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
        TestPattern(r"New AST:(.|\n)*?return A \* x",
                     "Test if optimization framework invokes SparseAccelerator to generate a new AST."
        ),
        exception_pattern
    ]
)

const spmv_sanity_test1 = Test(
    "spmv-sanity-test1",
    "spmv-sanity-test1.jl ../matrices/tiny-diag.mtx",
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

const CG_MATRIX = "../matrices/bcsstk14.mtx"

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
        TestPattern(r"New AST:(.|\n)*?mknobA.* = \(SparseAccelerator.new_matrix_knob\)\(:A,true,true,true,true,false,false\)",
                     "Test if mknobA is generated and is constant valued and constant structured"
        ),
        TestPattern(r"New AST:(.|\n)*?add_mknob_to_fknob\)\(.*mknobA.*,..*fknob.*\)",
                     "Test if mknobA is added to a function knob (for SpMV)"
        ),
        TestPattern(r"New AST:(.|\n)*?Ap = \(SparseAccelerator.SpMV\)\(1,A.*,p.*,.*fknob.*\)",
                     "Test if Ap = A * p has been replaced with SpMV with context info"
        ),
        TestPattern(r"New AST:(.|\n)*?set_reordering_decision_maker",
                     "Test reordering"
        ),
        TestPattern(r"\(SparseAccelerator.reordering\)\(__fknob\d*__,__reordering_status\d*__,A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobA\d*__,U,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobU\d*__,:__delimitor__,p,SparseAccelerator.ROW_PERM,r,SparseAccelerator.ROW_PERM,x,SparseAccelerator.ROW_PERM\)",
                     "Test reordering"
        ),      
        TestPattern(r"reverse_reordering\)\(__reordering_status\d*__,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
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
    string("context-test2-ilu0.jl  ../matrices/lap3d_4x4x4.mtx"),
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
        TestPattern(r"New AST:(.|\n)*?B = .*\(SparseAccelerator.ADB\)\(A.*,D.*,__AT__.*,__fknob\d*?__",
                     "Test if accelerated ipm-ref generates ADB"
        ),
        TestPattern(r"New AST:(.|\n)*?B = .*\(SparseAccelerator.ADB\).*\n.*propagate_matrix_info.*mknobB.*mknobExpr",
                     "Test if accelerated ipm-ref generates propagate_matrix_info after B = ADB"
        ),
        TestPattern(r"New AST:(.|\n)*?R = .*\(SparseAccelerator.cholfact_int32\)\(B.*,__fknob\d*?__",
                     "Test if accelerated ipm-ref generates cholfact_int32"
        ),
        TestPattern(r"New AST:(.|\n)*?dy = .*\(SparseAccelerator.cholfact_inverse_divide\)\(R.*,t2.*,__fknob\d*?__",
                     "Test if accelerated ipm-ref generates cholfact_inverse_divide"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.set_derivative\)\(.*mknobA.*,SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE,.*mknob__AT.*\)",
                     "Test if accelerated ipm-ref generates a derivative between A and A'"
        ),
        # TODO: enable these patterns
#        TestPattern(r"New AST:(.|\n)*?Rd = .*\(SparseAccelerator.SpMV\)\(1,__AT__,y,1,s-p,fknob_spmv.*\)",
#                     "Test call replacement. Currently cannot match as we need to handle AT specially."
#        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Rp,1,A,x,-1,b,0,.*?fknob.*?\)",
                     "Test call replacement."
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.element_wise_multiply!\)\(Rc,x,s\)",
                     "Test call replacement."
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.WAXPB!\)\(Rc,1,Rc,-\(\(Main.min\)\(0.1,100 \* mu.*?\).*? \* mu.*?\)",
                     "Test call replacement."
        ),
        #TestPattern(r"Original sum of x=715375.98850000",
                     #"Test original ipm-ref"
        #),
        TestPattern(r"Original iter 26, resid =  3.93e-15",
                     "Test original ipm-ref"
        ),
        #TestPattern(r"Accelerated sum of x=715375.98850000",
                     #"Test ipm-ref with context-sensitve optimization"
        #),
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
        TestPattern(r"Original iter 26, resid =  3.93e-15",
                     "Test original ipm-ref"
        ),
        TestPattern(r"Manual_context iter 26, resid =  3.93e-15",
                     "Test ipm-ref with context-sensitve optimization"
        ),
        exception_pattern
    ]
)

const context_test5 = Test(
    "context-test5",
    "context-test5.jl ../matrices/lap3d_4x4x4.mtx",
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
        TestPattern(r"New AST:(.|\n)*?mknobA.* = \(SparseAccelerator.new_matrix_knob\)\(:A,true,true,true,true,false,false\)",
                     "Test if mknobA is generated and is constant valued, constant structured, symmetric valued and structured"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobL.* = \(SparseAccelerator.new_matrix_knob\)\(:L,true,true,false,false,false,false\)",
                     "Test if mknobL is generated and is constant valued and constant structured"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobU.* = \(SparseAccelerator.new_matrix_knob\)\(:U,true,true,false,false,false,false\)",
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
#        TestPattern(r"New AST:(.|\n)*?Ap = \(SparseAccelerator.SpMV\)\(1,A,p,.*fknob.*\)",
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Ap,A,p,__fknob\d*?__\)",
                     "Test if Ap = A * p has been replaced with SpMV with context info"
        ),
        TestPattern(r"New AST:(.|\n)*?set_reordering_decision_maker",
                     "Test reordering"
        ),
#        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.reordering\)\(__fknob\d*?__,##reordering_status#\d*?,(U,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,|A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,){2,2}:__delimitor__,(r,SparseAccelerator.ROW_PERM,?|p,SparseAccelerator.ROW_PERM,?|x,SparseAccelerator.ROW_PERM,?){3,3}\)",
        TestPattern( r"New AST:(.|\n)*?SparseAccelerator.reordering\)\(__fknob\d*?__,__reordering_status\d*?__,(U,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobU\d*?__,|A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobA\d*?__,){2,2}:__delimitor__,(r,SparseAccelerator.ROW_PERM,?|p,SparseAccelerator.ROW_PERM,?|x,SparseAccelerator.ROW_PERM,?|Ap,SparseAccelerator.ROW_PERM,?){4,4}\)",
                      "Test reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?reverse_reordering\)\(__reordering_status\d*?__,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
                     "Test reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?r = b - \(SparseAccelerator.SpMV\)\(A,x\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?normr0 = .*SparseAccelerator.norm\)\(r.*\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?rz = .*SparseAccelerator.dot\)\(r.*,z.*\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?alpha = old_rz.*?/.*?SparseAccelerator.dot\)\(p.*?,Ap.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.WAXPBY!\)\(x,1,x.*?,alpha.*?,p.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.WAXPBY!\)\(r,1,r.*?,-alpha.*?,Ap.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?rel_err = .*?SparseAccelerator.norm\)\(r.*?\).*?/ normr0",
                     "Test call replacement"
        ),
# TODO: enable this test.
#        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:copy!\)\(z.*?,r.*?\)",
#                     "Test call replacement"
#        )
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.fwdTriSolve!\)\(L,z,.*?fknob.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.bwdTriSolve!\)\(U,z,.*?fknob.*?\)",
                     "Test call replacement"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.WAXPBY!\)\(p,1,z.*?,beta.*?,p.*?\)",
                     "Test call replacement"
        ),
        exception_pattern
    ]
)

const pagerank_test1 = Test(
    "pagerank-test1",
    "pagerank.jl  ../matrices/hmatrix.1024.mtx",
    [
        TestPattern(r"Original:(.|\n)*?\s*error = 1.53\d*e-8",
                     "Test original pagerank"
        ),
        TestPattern(r"Accelerated:(.|\n)*?\s*error = 1.53\d*e-8",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?set_reordering_decision_maker",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Ap,1 - r.*?,A.*?,p.*?,0,p.*?,r.*?,d_inv.*?,__fknob\d*?__\)",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"reverse_reordering\)\(__reordering_status\d*?__,:__delimitor__,p,SparseAccelerator.ROW_PERM\)",
                     "Test pagerank with reordering"
        ),
        TestPattern(r"err = .*?\(SparseAccelerator.norm\)\(Ap.*? - p.*?\).*? / .*?\(SparseAccelerator.norm\)\(p.*?\)",
                     "Test pagerank if call replacement of norm happens"
        ),
        exception_pattern
    ]
)

const cosp2_test1 = Test(
    "cosp2-test1",
    "cosp2.jl ../matrices/hmatrix.512.mtx",
    [
        TestPattern(r"Original:\nX sum = 6106.4049\d*, max = 1.2808\d*\nNumber of iterations = 25(.|\n)*?End original.",
                     "Test original"
        ),
        TestPattern(r"CoSP2_call_replacement:\nX sum = 6106.4049\d*, max = 1.2808\d*\nNumber of iterations = 25(.|\n)*?End CoSP2_call_replacement.",
                     "Test CoSP2_call_replacement"
        ),
        TestPattern(r"CoSP2_call_replacement_and_context_opt:\nX sum = 6106.4049\d*, max = 1.2808\d*\nNumber of iterations = 25(.|\n)*?End CoSP2_call_replacement_and_context_opt.",
                     "Test CoSP2_call_replacement_and_context_opt"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobX.*? = \(SparseAccelerator.new_matrix_knob\)\(:X,false,.*?,true,true,false,false\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.add_mknob_to_fknob\)\(.*?mknobX.*?,.*?fknob.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?trX = \(SparseAccelerator.trace\)\(X.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?X2 = \(SparseAccelerator.SpSquareWithEps\)\(X.*?,eps.*?,.*?fknob.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?trX2 = \(SparseAccelerator.trace\)\(X2.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?X = \(SparseAccelerator.SpAdd\)\(2,X.*?,-1,X2.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?X sum = 6106.4049\d*, max = 1.2808\d*\nNumber of iterations = 25(.|\n)*?End accelerated",
                     "Test accelerated"
        ),
        exception_pattern
    ]
)

#const MTX_covtype = "../matrices/covtype.mtx" 

if false
const lbfgs_test1 = Test(
    "lbfgs-test1",
    "lbfgs.jl $MTX_covtype",
    [
        TestPattern(r"Original L-BFGS:      32 iterations f = 0.00000000004134",
                     "Test original"
        ),
        TestPattern(r"Opt L-BFGS:      32 iterations f = 0.00000000004134",
                     "Test optimized version"
        ),
        TestPattern(r"Opt_with_reordering L-BFGS:      33 iterations f = 0.0000000000",
                     "Test Opt_with_reordering"
        ),
        TestPattern(r"Accelerated L-BFGS:      3[345] iterations f = 0.0000000",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobXt.*? = \(SparseAccelerator.new_matrix_knob\)\(:Xt,true,true,false,false,false,false\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobX.*? = \(SparseAccelerator.new_matrix_knob\)\(:X,true,true,false,false,false,false\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.set_derivative\)\(__mknobXt\d*?__,SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE,__mknobX\d*?__\)",
                     "Test accelerated"
        ),                                 
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Xw,X,x,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.element_wise_multiply!\)\(yXw,y,Xw\)",
                     "Test accelerated"
        ),           
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.abs!\)\(__temp.*?,yXw\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(__temp.*?,-1,__temp.*?,0,__temp.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.exp!\)\(__temp.*?,__temp.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.log1p!\)\(__temp.*?,__temp.*?\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.sum\)\(\(SparseAccelerator.WAXPBY!\)\(__temp.*?,1,__temp.*?,-1,\(SparseAccelerator.min!\)\(__temp.*?,yXw,0\)\)\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?fk0 = s / m \+ \(lambda / 2\) \* \(SparseAccelerator.dot\)\(x,x\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.exp!\)\(__temp.*?,yXw\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.element_wise_divide!\)\(temp,y,\(SparseAccelerator.WAXPB!\)\(__temp.*?,1,__temp.*?,1\)\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(dfk,-1 / m,Xt,temp,lambda,x,0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.reordering\)\(__fknob\d*?__,__reordering_status\d*?__,X,SparseAccelerator.COL_INV_PERM,SparseAccelerator.ROW_INV_PERM,__mknobX\d*?__,:__delimitor__,y,SparseAccelerator.COL_INV_PERM,dx,SparseAccelerator.ROW_PERM\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?unless \(SparseAccelerator.norm\)\(dfk\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.reverse_reordering\)\(__reordering_status\d*?__,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(w,1,x,-alpha,dfk\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Xw,1,X,w,0.0,Xw,0.0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.element_wise_multiply!\)\(yXw,y,Xw\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?fk = s / m \+ \(lambda / 2\) \* \(SparseAccelerator.dot\)\(w,w\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.dot\)\(dfk,dfk\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(x,1,x,alpha,dx\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Xw,1,X,x,0.0,Xw,0.0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(dfkp1,-1 / m,Xt,temp,lambda,x,0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(__temp.*?,1,dfkp1,-1,dfk\)",
                     "Test accelerated"
        ),    
        exception_pattern
    ]
)

const lbfgs_test2 = Test(
    "lbfgs-test2",
    "lbfgs-new.jl $MTX_covtype",
    [
        TestPattern(r"Original L-BFGS:      32 iterations f = 0.00000000004134",
                     "Test original"
        ),
        TestPattern(r"Opt L-BFGS:      32 iterations f = 0.00000000004134",
                     "Test optimized version"
        ),
        TestPattern(r"Opt_with_reordering L-BFGS:      33 iterations f = 0.0000000000",
                     "Test Opt_with_reordering"
        ),
        TestPattern(r"Accelerated L-BFGS:      3[345] iterations f = 0.00000000",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobXt.*? = \(SparseAccelerator.new_matrix_knob\)\(:Xt,true,true,false,false,false,false\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?mknobX.*? = \(SparseAccelerator.new_matrix_knob\)\(:X,true,true,false,false,false,false\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.set_derivative\)\(__mknobXt\d*?__,SparseAccelerator.DERIVATIVE_TYPE_TRANSPOSE,__mknobX\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Xw,1,X,x,0.0,Xw,0.0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.element_wise_multiply!\)\(yXw,y,Xw\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(dfk,-1 / m,Xt,temp,lambda,x,0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.reordering\)\(__fknob\d*?__,__reordering_status\d*?__,X,SparseAccelerator.COL_INV_PERM,SparseAccelerator.ROW_INV_PERM,__mknobX\d*?__,:__delimitor__,y,SparseAccelerator.COL_INV_PERM,dx,SparseAccelerator.ROW_PERM\)
",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?unless \(SparseAccelerator.norm\)\(dfk\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.reverse_reordering\)\(__reordering_status\d*?__,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(w,1,x,-alpha,dfk\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(Xw,1,X,w,0.0,Xw,0.0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.element_wise_multiply!\)\(yXw,y,Xw\)",
                     "Test accelerated"
        ),        
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator.dot\)\(dfk,dfk\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(dfkp1,-1 / m,Xt,temp,lambda,x,0,__fknob\d*?__\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.reverse_reordering\)\(__reordering_status\d*?__,:__delimitor__,x,SparseAccelerator.ROW_PERM\)",
                     "Test accelerated"
        ),
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(__temp.*?,1,dfkp1,-1,dfk\)",
                     "Test accelerated"
        ),    
        exception_pattern
    ]
)
end

const MTX_small_diag = "../matrices/small-diag.mtx"

const liveness_test1 = Test(
    "liveness-test1",
    "liveness-test1.jl $MTX_small_diag",
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
    "call-replacement-test1.jl $MTX_small_diag",
    [
        TestPattern(r"AST:(.|\n)*?Main.dot(.|\n)*?New AST(.|\n)*?SparseAccelerator.dot",
                     "Test call replacement of Main.dot with SparseAccelerator.dot."
        ),
        exception_pattern
    ]
)

const call_replacement_test2 = Test(
    "call-replacement-test2",
    "call-replacement-test2.jl $MTX_small_diag",
    [
        TestPattern(r"AST:(.|\n)*?\(SparseAccelerator.SpMV\)\(A,x\)",
                     "Test call replacement of * with SparseAccelerator.SpMV."
        ),
        exception_pattern
    ]
)

const call_replacement_test3 = Test(
    "call-replacement-test3",
    "call-replacement-test3.jl $MTX_small_diag",
    [
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(y,A,x\)",
                     "Test call replacement of SpMV! for A_mul_B!(y, A, x)."
        ),
        exception_pattern
    ]
)

const call_replacement_test4 = Test(
    "call-replacement-test4",
    "call-replacement-test4.jl $MTX_small_diag",
    [
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.SpMV!\)\(y,0.1,A,x,0.1,y,0.0\)",
                     "Test call replacement of SpMV for A_mul_B!(0.1, A, x, 0.1, y)."
        ),
        exception_pattern
    ]
)

const call_replacement_test5 = Test(
    "call-replacement-test5",
    "call-replacement-test5.jl $MTX_small_diag",
    [
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(x,1,x,alpha,p\)",
                     "Test call replacement of WAXPBY! for x += alpha * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test6 = Test(
    "call-replacement-test6",
    "call-replacement-test6.jl $MTX_small_diag",
    [
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(x,1,x,-alpha,p\)",
                     "Test call replacement of WAXPBY! for x -= alpha * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test7 = Test(
    "call-replacement-test7",
    "call-replacement-test7.jl $MTX_small_diag",
    [
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(p,1,r,beta,p\)",
                     "Test call replacement of WAXPBY! for p = r + beta * p."
        ),
        exception_pattern
    ]
)

const call_replacement_test8 = Test(
    "call-replacement-test8",
    "call-replacement-test8.jl $MTX_small_diag",
    [
        TestPattern(r"New AST:(.|\n)*?SparseAccelerator,:SpMV\!\)\(p,1 - r::Float64::Float64,A::[^\s\\]*SparseMatrixCSC\{Float64,Int32\},p::Array\{Float64,1\},0,p::Array\{Float64,1\},r::Float64\)",
                     "Test call replacement of SpMV! in simple page rank."
        ),
        exception_pattern
    ]
)

const call_replacement_test9 = Test(
    "call-replacement-test9",
    "call-replacement-test9.jl $MTX_small_diag",
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
    "call-replacement-test10.jl $MTX_small_diag",
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
    "call-replacement-test11.jl $MTX_small_diag",
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
        TestPattern(r"New AST:(.|\n)*?\(SparseAccelerator.WAXPBY!\)\(p,1,z,beta,p\)",
                     "Test call replacement of WAXPBY! for p = r + beta * p."
        ),
        exception_pattern
    ]
)

const name_resolution_test1 = Test(
    "name-resolution-test1",
    "name-resolution-test1.jl $MTX_small_diag",
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

const single_def_test2 = Test(
    "single-def-test2",
    "single-def-test2.jl",
    [
        TestPattern(Regex("Single-defs discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test ipm-ref that A B are recognized as single-defs in the loop."
        ),
        exception_pattern
    ]
)


const set_matrix_property_test1 = Test(
    "set-matrix-property-test1",
    "set-matrix-property-test1.jl",
    [
        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that R is recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const set_matrix_property_test2 = Test(
    "set-matrix-property-test2",
    "set-matrix-property-test2.jl",
    [

        TestPattern(Regex("Loop0 Matrix structures discovered:.*\\n.*" * gen_set_regex_string([:A, :L, :U])),
                     "Test that A L U are recognized as matrics in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :L, :U])),
                     "Test that A L U are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that R is recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const set_matrix_property_test3 = Test(
    "set-matrix-property-test3",
    "set-matrix-property-test3.jl",
    [
        TestPattern(Regex("Func Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Func Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Func Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that R is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that A B are recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as constant in structure."
        ),

        TestPattern(Regex("Loop0 Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that R is recognized as constant in structure."
        ),
        exception_pattern
    ]
)


const set_matrix_property_test4 = Test(
    "set-matrix-property-test4",
    "set-matrix-property-test4.jl",
    [

        TestPattern(Regex("Func Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as matrics in structure."
        ),

        TestPattern(Regex("Loop0 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as matrics in structure."
        ),

        TestPattern(Regex("Loop3 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that A B are recognized as constant in structure."
        ),

        exception_pattern
    ]
)

const set_matrix_property_test5 = Test(
    "set-matrix-property-test5",
    "set-matrix-property-test5.jl",
    [
        TestPattern(r"Func Upper/Lower matrix discovered:.*\n[\s]*A is lower of B",
                     "Test that A is recognized as lower part of B."
        ),
        TestPattern(r"Loop0 Upper/Lower matrix discovered:.*\n[\s]*A is lower of B",
                     "Test that A is recognized as lower part of B."
        ),
        TestPattern(r"Loop3 Upper/Lower matrix discovered:.*\n[\s]*A is lower of B.*\n[\s]*B is upper of C",
                     "Test that A is recognized as matrics in structure."
        ),
        exception_pattern
    ]
)

const set_matrix_property_test6 = Test(
    "set-matrix-property-test6",
    "set-matrix-property-test6.jl",
    [

        TestPattern(Regex("Func Structure only discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as structure-only matrix."
        ),

        TestPattern(Regex("Loop0 Structure only discovered:.*\\n.*" * gen_set_regex_string([:A, :C])),
                     "Test that A C are recognized as structure-only matrics."
        ),

        exception_pattern
    ]
)

const set_matrix_property_test7 = Test(
    "set-matrix-property-test7",
    "set-matrix-property-test7.jl",
    [
        TestPattern(r"Func Transpose matrix discovered:.*\n[\s]*A is transpose of B",
                     "Test that A is recognized as transpose of B."
        ),
        TestPattern(r"Loop0 Transpose matrix discovered:.*\n[\s]*A is transpose of B\n[\s]*B is transpose of C",
                     "Test that A is recognized as transpsoe of C."
        ),

        exception_pattern
    ]
)


const prop_constant_structure_test1 = Test(
    "prop-constant-structure-test1",
    "prop-constant-structure-test1.jl ipm/mps/osa-14",
    [
        TestPattern(Regex("Loop-4 Constant structures discovered:.*\\n.*" * gen_set_regex_string([:A, :B, :D, :R])),
                     "Test that A B D R are recognized as constant in structure."
        ),
        exception_pattern
    ]
)

const prop_symmetric_value_test1 = Test(
    "prop-symmetric-value-test1",
    "prop-symmetric-value-test1.jl",
    [
        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:B, :E, :F, :G, :S])),
                     "Test that B E F S is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const prop_symmetric_value_test2 = Test(
    "prop-symmetric-value-test2",
    "prop-symmetric-value-test2.jl",
    [
        TestPattern(r"Value symmetry discovered:.*\n.*\[.*:A.*\]",
                     "Test that A is recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const prop_symmetric_value_test3 = Test(
    "prop-symmetric-value-test3",
    "prop-symmetric-value-test3.jl",
    [
        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:B, :C, :G, :S])),
                     "Test that B E F G S are recognized as symmetric in value."
        ),
        exception_pattern
    ]
)

const prop_symmetric_value_test4 = Test(
    "prop-symmetric-value-test4",
    "prop-symmetric-value-test4.jl",
    [
        TestPattern(Regex("Loop0 Value symmetry discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that B E F G S are recognized as symmetric in value."
        ),
        exception_pattern
    ]
)


const prop_symmetric_structure_test1 = Test(
    "prop-symmetric-structure-test1",
    "prop-symmetric-structure-test1.jl",
    [
        TestPattern(Regex("Structure symmetry discovered:.*\\n.*" * gen_set_regex_string([:A, :B])),
                     "Test that A B are recognized as symmetric in structure."
        ),
        exception_pattern
    ]
)

const prop_structure_only_test1 = Test(
    "prop-structure-only-test1",
    "prop-structure-only-test1.jl",
    [
        TestPattern(Regex("Func Structure only discovered:.*\\n.*" * gen_set_regex_string([:B, :C, :E])),
                     "Test that B C E are recognized as pattern-only structure."
        ),
        exception_pattern
    ]
)


const prop_structure_only_test2 = Test(
    "prop-structure-only-test2",
    "prop-structure-only-test2.jl",
    [
        TestPattern(Regex("Func Structure only discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as pattern-only structure."
        ),
        exception_pattern
    ]
)

const prop_structure_only_test3 = Test(
    "prop-structure-only-test3",
    "prop-structure-only-test3.jl",
    [
        TestPattern(Regex("Func Structure only discovered:.*\\n.*" * gen_set_regex_string([:A])),
                     "Test that A is recognized as pattern-only structure."
        ),
        exception_pattern
    ]
)


const prop_lower_upper_test1 = Test(
    "prop-lower-upper-test1",
    "prop-lower-upper-test1.jl",
    [
        TestPattern(r"Func Upper/Lower matrix discovered:.*\n[\s]*D is upper of A.*\n[\s]*GenSym\(0\) is upper of A",
                     "Test D and GenSym(0) are upper part of A."
        ),
        exception_pattern
    ]
)


const prop_lower_upper_test2 = Test(
    "prop-lower-upper-test2",
    "prop-lower-upper-test2.jl",
    [
        TestPattern(r"Loop0 Upper/Lower matrix discovered:.*\n[\s]*C is lower of A.*\n[\s]*D is upper of A.*\n[\s]*E is lower of A.*\n[\s]*F is upper of A.*\n[\s]*L is lower of A",
                     "Test that A B are recognized as symmetric in structure."
        ),
        exception_pattern
    ]
)


const prop_transpose_test1 = Test(
    "prop-transpose-test1",
    "prop-transpose-test1.jl",
    [
        TestPattern(r"Func Transpose matrix discovered:.*\n[\s]*At is transpose of A.*\n[\s]*Att is transpose of At.*\n[\s]*B is transpose of A.*",
                     "Test D and GenSym(0) are upper part of A."
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
#    lbfgs_test1,
#    lbfgs_test2,
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
    single_def_test2,
    constant_value_test1,
    set_matrix_property_test1,
    set_matrix_property_test2,
    set_matrix_property_test3,
    set_matrix_property_test4,
    set_matrix_property_test5,
    set_matrix_property_test6,
    set_matrix_property_test7,
    prop_constant_structure_test1,
    prop_symmetric_value_test1,
    prop_symmetric_value_test2,
    prop_symmetric_value_test3,
    prop_symmetric_value_test4,
#    symmetric_structure_test1,
    prop_structure_only_test1,
    prop_structure_only_test2,
    prop_structure_only_test3,
    prop_lower_upper_test1,
    prop_lower_upper_test2,
    prop_transpose_test1,
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
#    lbfgs_test1
]

# If true, use pcregrep for regular expression match. 
# If false, use Julia (If PCRE JIT stack overflow: enlarge JIT_STACK_MAX_SIZE in
# julia/base/pcre.jl (times it with 10) and retry).
const USE_PCREGREP_REGEX_MATCH = true

function get_julia_ver()
    s, p = open(`$julia_command -v`)
    readline(s)
end

if !isreadable(julia_command)
    error("Please install (softlink) julia command to \"" * julia_command  * "\".")
elseif !ismatch(r"\.*0.4.1", get_julia_ver())
#    error("Wrong julia version! 0.4.1 is required!")
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
passed = Dict()
if nprocs() == 1 # single thread mode
    for test in tests
        print("Testing ", test.name)
        s = run_test(test)
        if s
            succ = succ + 1
            println(": Pass")
	    passed[test.name] = true
        else
            println(": FAIL. See ", test.name * ".log")
	    passed[test.name] = false
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

end

# Print out the results in a fixed order for easier checking
println("\nFinal report:")
for test in tests
    println(test.name, ": ", passed[test.name] ? "Pass" : "FAIL. See $(test.name).log")
end

println("Total: ", total)
println("Pass : ", succ)
println("Fail : ", total-succ)

