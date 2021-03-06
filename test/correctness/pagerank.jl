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

include("../../src/Sparso.jl")
include("../../src/simple-show.jl")
using Sparso

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)

function pagerank(A, p, r, d_inv, maxiter) # p: initial rank, r: damping factor
#  set_matrix_property(Dict(
#    :A => SA_SYMM_STRUCTURED | SA_SYMM_VALUED | SA_STRUCTURE_ONLY))

  bytes = maxiter*(nnz(A)*4 + size(A, 1)*4*8)

  p = p.*d_inv
  d_inv = copy(d_inv)
  Ap = zeros(size(A, 1))

  t = time()

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

function pagerank_with_init(A0, r, d_inv, maxiter) # p: initial rank, r: damping factor
  set_matrix_property(Dict(
    :A0 => SA_SYMM_VALUED)
  )

  A0 = spones(A0)

  m = size(A0, 1)
  p = repmat([1/m], m)

  d = max(convert(Array{eltype(A0),1}, vec(sum(A0, 2))), 1) # num of neighbors
  A = scale(A0,1./d)

  bytes = maxiter*(nnz(A)*4 + m*4*8)

  p = p.*d_inv
  d_inv = copy(d_inv)
  Ap = zeros(size(A, 1))

  t = time()

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


function pagerank_call_replacement(A, p, r, d_inv, maxiter, do_print) # p: initial rank, r: damping factor
  bytes = maxiter*(nnz(A)*4 + size(A, 1)*4*8)
  A = copy(A)
  p = p.*d_inv
  d_inv = copy(d_inv)
  Ap = copy(p)

  t = time()
  reorder_time = 0.

  __mknobA = (Sparso.new_matrix_knob)(:A, true, true, true, true, false, false)

  fknob_spmv = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobA, fknob_spmv)

  for i = 1:maxiter
    #Ap = ((1-r) *A * p + r).*d_inv
    Sparso.SpMV!(Ap, 1 - r, A, p, 0, p, r, d_inv, fknob_spmv)
    if do_print && i == maxiter
      err = norm(Ap - p)/norm(p)
      println("error = $err")
    end

    temp = Ap
    Ap = p
    p = temp
  end

  (Sparso.delete_matrix_knob)(__mknobA)
  (Sparso.delete_function_knob)(fknob_spmv)

  t = time() - t
  if do_print
    println("pagerank takes $t sec ($(bytes/t/1e9) gbps)")
  end

  p
end

function pagerank_call_replacement_and_context_opt(A, p, r, d_inv, maxiter, do_print) # p: initial rank, r: damping factor
  bytes = maxiter*(nnz(A)*4 + size(A, 1)*4*8)
  A = copy(A)
  p = p.*d_inv
  d_inv = copy(d_inv)
  Ap = copy(p)

  t = time()
  reorder_time = 0.

  __mknobA = (Sparso.new_matrix_knob)(:A, true, true, true, true, true, false)

  fknob_spmv = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobA, fknob_spmv)

  for i = 1:maxiter
    #Ap = ((1-r) *A * p + r).*d_inv
    Sparso.SpMV!(Ap, 1 - r, A, p, 0, p, r, d_inv, fknob_spmv)
    if do_print && i == maxiter
      err = norm(Ap - p)/norm(p)
      println("error = $err")
    end

    temp = Ap
    Ap = p
    p = temp
  end

  (Sparso.delete_matrix_knob)(__mknobA)
  (Sparso.delete_function_knob)(fknob_spmv)

  t = time() - t
  if do_print
    println("pagerank takes $t sec ($(bytes/t/1e9) gbps)")
  end

  p
end

function pagerank_call_replacement_and_context_opt_and_reordering(A, p, r, d_inv, maxiter, do_print) # p: initial rank, r: damping factor
  bytes = maxiter*(nnz(A)*4 + size(A, 1)*4*8)
  A = copy(A)
  p = p.*d_inv
  d_inv = copy(d_inv)
  Ap = copy(p)

  t = time()
  reorder_time = 0.

  __mknobA = (Sparso.new_matrix_knob)(:A, true, true, true, true, true, false)

  fknob_spmv = (Sparso.new_function_knob)()
  (Sparso.add_mknob_to_fknob)(__mknobA, fknob_spmv)

  (Sparso.set_reordering_decision_maker)(fknob_spmv)
  reordering_status = [false, C_NULL, C_NULL, C_NULL, C_NULL, reorder_time]

  for i = 1:maxiter
    #Ap = ((1-r) *A * p + r).*d_inv
    Sparso.SpMV!(Ap, 1 - r, A, p, 0, p, r, d_inv, fknob_spmv)
    if do_print && i == maxiter
      err = norm(Ap - p)/norm(p)
      println("error = $err")
    end

    temp = Ap
    Ap = p
    p = temp
  end

  (Sparso.delete_matrix_knob)(__mknobA)
  (Sparso.delete_function_knob)(fknob_spmv)
  (Sparso.reverse_reordering)(reordering_status, :__delimitor__, p, Sparso.ROW_PERM)

  t = time() - t
  if do_print
    println("pagerank takes $t sec ($(bytes/t/1e9) gbps)")
  end

  p
end

A0 = matrix_market_read(ARGS[1], true, true)
A0 = spones(A0)

m = size(A0, 1)
p = repmat([1/m], m)
r = 0.15

d = max(convert(Array{eltype(A0),1}, vec(sum(A0, 2))), 1) # num of neighbors
A = scale(A0,1./d)

maxiter = 100
bytes = maxiter*(nnz(A)*4 + m*4*8)

d_inv = 1./d
x = pagerank(A0, p, r, d_inv, maxiter)
println("\nOriginal: ")
x = pagerank(A0, p, r, d_inv, maxiter)
println("End original.")

println("\nManual_call_replacement:")
x = pagerank_call_replacement(A0, p, r, d_inv, maxiter, false)

Sparso.reset_spmp_spmv_time()
Sparso.reset_knob_spmv_time()
x = pagerank_call_replacement(A0, p, r, d_inv, maxiter, true)
t = Sparso.get_spmp_spmv_time()
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = Sparso.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")
println("End Manual_call_replacement.")

println("\nManual_call_replacement_and_context_opt:")
x = pagerank_call_replacement_and_context_opt(A, p, r, d_inv, maxiter, false)

Sparso.reset_spmp_spmv_time()
Sparso.reset_knob_spmv_time()
#Sparso.set_knob_log_level(1)
x = pagerank_call_replacement_and_context_opt(A0, p, r, d_inv, maxiter, true)
t = Sparso.get_spmp_spmv_time()
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = Sparso.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")

println("End Manual_call_replacement_and_context_opt.")

println("\nManual_call_replacement_and_context_opt_and_reorder:")
x = pagerank_call_replacement_and_context_opt_and_reordering(A, p, r, d_inv, maxiter, false)

Sparso.reset_spmp_spmv_time()
Sparso.reset_knob_spmv_time()
#Sparso.set_knob_log_level(1)
x = pagerank_call_replacement_and_context_opt_and_reordering(A0, p, r, d_inv, maxiter, true)
t = Sparso.get_spmp_spmv_time()
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = Sparso.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")

println("End Manual_call_replacement_and_context_opt_and_reorder.")

# copy A since we change A in-place
#Sparso.set_knob_log_level(0)
A2 = copy(A0)
#@acc x= pagerank(A2, p, r, d_inv, maxiter)
@acc x= pagerank_with_init(A2, r, d_inv, maxiter)

println("\nAccelerated: ")
Sparso.reset_spmp_spmv_time()
Sparso.reset_knob_spmv_time()
#Sparso.set_knob_log_level(1)
# It seems that p has been changed by sparse accelerator
p = repmat([1/m], m)
#@acc x = pagerank(A0, p, r, d_inv, maxiter)
@acc x = pagerank_with_init(A2, r, d_inv, maxiter)
t = Sparso.get_spmp_spmv_time()
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = Sparso.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")
println("End accelerated.")
