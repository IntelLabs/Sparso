include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)

function pagerank(A, p, r, maxiter) # p: initial rank, r: damping factor
  set_matrix_property(Dict(
    :A => SA_SYMM_STRUCTURED | SA_SYMM_VALUED))

  bytes = maxiter*(nnz(A)*12 + size(A, 1)*3*8)

  t = time()

  Ap = zeros(size(A, 1))

  for i = 1:maxiter
    A_mul_B!(1 - r, A, p, 0, Ap)
    Ap = Ap + r
    #Ap = (1-r)*A*p + r

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
A = scale(A,1./d)

maxiter = 100

p2 = copy(p)
x = pagerank(A, p2, r, maxiter)
println("Original: ")
p2 = copy(p)
x = pagerank(A, p2, r, maxiter)

# copy A since we change A in-place
A2 = copy(A)
p2 = copy(p)
@acc x= pagerank(A2, p2, r, maxiter)

println("\nAccelerated: ")
SparseAccelerator.reset_spmp_spmv_time()
SparseAccelerator.reset_knob_spmv_time()
@acc x= pagerank(A, p, r, maxiter)
t = SparseAccelerator.get_spmp_spmv_time()
bytes = maxiter*(nnz(A)*12 + m*3*8)
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = SparseAccelerator.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")
