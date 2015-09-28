include("../../src/SparseAccelerator.jl")
using SparseAccelerator

maxiter = 100

function pagerank(A, p, r) # p: initial rank, r: damping factor
  bytes = maxiter*(nnz(A)*12 + size(A, 1)*3*8)
  p = copy(p)
  Ap = copy(p)

  t = time()
  reorder_time = 0.

  __mknobA = (SparseAccelerator.new_matrix_knob)(A, true, true, true, true, false, false)

  fknob_spmv = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(__mknobA, fknob_spmv)

  (SparseAccelerator.set_reordering_decision_maker)(fknob_spmv)
  reordering_status = [false, C_NULL, C_NULL, C_NULL, C_NULL, reorder_time]

  for i = 1:maxiter
    #Ap = (1-r) *A * p + r
    SparseAccelerator.SpMV!(Ap, 1 - r, A, p, 0, p, r, fknob_spmv)
    SparseAccelerator.reordering(fknob_spmv, reordering_status, :__delimitor__)
    if i == maxiter
      err = norm(Ap - p)/norm(p)
      println("error = $err")
    end

    temp = Ap
    Ap = p
    p = temp
  end

  (SparseAccelerator.delete_matrix_knob)(__mknobA)
  (SparseAccelerator.delete_function_knob)(fknob_spmv)
  (SparseAccelerator.reverse_reordering)(reordering_status, :__delimitor__, p, SparseAccelerator.ROW_PERM)

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

x = pagerank(A, p, r)

println("\nManually optimized: ")
SparseAccelerator.reset_spmp_spmv_time()
SparseAccelerator.reset_knob_spmv_time()
x= pagerank(A, p, r)
t = SparseAccelerator.get_spmp_spmv_time()
bytes = maxiter*(nnz(A)*12 + m*3*8)
println("time spent on spmp spmv $t sec ($(bytes/t/1e9) gbps)")
t = SparseAccelerator.get_knob_spmv_time()
println("time spent on knob spmv $t sec ($(bytes/t/1e9) gbps)")
