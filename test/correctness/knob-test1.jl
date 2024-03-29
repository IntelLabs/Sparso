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

# $ julia knob-test1.jl ../matrices/bcsstk14.mtx 
# assymmetric (2, 8) exists but (8, 2) doesn't
# julia: SymGS.cpp:305: bool SpMP::getSymmetricNnzPattern(const SpMP::CSR *, int **, int **, int **, int **): Assertion `sym.isSymmetric(false, true)' failed.
# Aborted (core dumped)

include("../../src/Sparso.jl")
include("../../src/simple-show.jl")
include("./utils.jl")
using Sparso

function pcg_symgs(x, A, b)
    L, U = Sparso.ilu(A)
    r = b - Sparso.SpMV(1,A,x)
    z = copy(r)
    Sparso.fwdTriSolve!(z,L,r, C_NULL)
end

A = matrix_market_read(ARGS[1], true, true)
m = size(A, 1)
b = ones(Float64, m)
x = zeros(Float64, m)
pcg_symgs(x, A, b)
