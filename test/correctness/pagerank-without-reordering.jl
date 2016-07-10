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

include("../../src/SparseAccelerator.jl")
include("../../src/simple-show.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS)

function pagerank(A, p, r, maxiter) # p: initial rank, r: damping factor
  set_matrix_property(Dict(
    :A => SA_SYMM_STRUCTURED | SA_SYMM_VALUED | SA_STRUCTURE_ONLY))

  bytes = maxiter*(nnz(A)*4 + size(A, 1)*4*8)

  p = p.*d_inv
  d_inv = copy(d_inv)
  Ap = zeros(size(A, 1))

  t = time()

  Ap = zeros(size(A, 1))

  for i = 1:maxiter
    Ap = (1-r)*A*p + r
    Ap = Ap.*d_inv

    if i == maxiter
      err = norm(Ap - p)/norm(p)
      println("error = $err")
    end

    temp = Ap
    Ap = p
    p = temp
  end

  t = time() - t
  println("pagerank takes $t sec ($(bytes/t/1e9) gbps)")

  p
end

A = matrix_market_read(ARGS[1], true, true)
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
r = 0.15

d = max(convert(Array{eltype(A),1}, vec(sum(A, 2))), 1) # num of neighbors
d_inv = 1./d

maxiter = 100

p2 = copy(p)
x = pagerank(A, p2, r, maxiter)
println("Original: ")
p2 = copy(p)
x = pagerank(A, p2, r, maxiter)

p2 = copy(p)
@acc x= pagerank(A, p2, r, maxiter)

println("\nAccelerated without reordering: ")
SparseAccelerator.reset_spmp_spmv_time()
SparseAccelerator.reset_knob_spmv_time()
@acc x= pagerank(A, p, r, maxiter)
t = SparseAccelerator.get_spmp_spmv_time()
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = SparseAccelerator.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")
