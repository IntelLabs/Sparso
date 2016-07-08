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

function generate_symmetric_sparse_matrix(m)
    A = SparseMatrixCSC{Cdouble, Cint}(sprand(m, m, 0.1))
    A = A*A'
    return A
end

function generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    A = generate_symmetric_sparse_matrix(m)
    
    # Make digonal nonzero
    for i = 1:m 
        if abs(A[i, i]) < 0.5
            A[i, i] = 1.0
        end 
    end
    
    return A
end

function check_symmetry(A)
    println("**** Checking symmetry")
    m = size(A, 1)
    for i = 1:m 
        for j = 1:m 
            if A[i, j] != A[j, i]
                println("Matrix is asymmetric!")
                println("A[", i, ",", j, "] != A[", j, ",", i, "]")
                break
            end
        end
    end
    println("Done checking")
    flush(STDOUT::IO)
end
