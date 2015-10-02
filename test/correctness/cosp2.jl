include("../../src/SparseAccelerator.jl")
using SparseAccelerator

function spmatmul_witheps{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}, eps;
                         sortindices::Symbol = :sortcols)
    mA, nA = size(A)
    mB, nB = size(B)
    nA==mB || throw(DimensionMismatch())

    colptrA = A.colptr; rowvalA = A.rowval; nzvalA = A.nzval
    colptrB = B.colptr; rowvalB = B.rowval; nzvalB = B.nzval
    # TODO: Need better estimation of result space
    nnzC = min(mA*nB, length(nzvalA) + length(nzvalB))
    colptrC = Array(Ti, nB+1)
    rowvalC = Array(Ti, nnzC)
    nzvalC = Array(Tv, nnzC)

    @inbounds begin
        ip = 1
        xb = zeros(Ti, mA)
        x  = zeros(Tv, mA)
        for i in 1:nB
            if ip + mA - 1 > nnzC
                resize!(rowvalC, nnzC + max(nnzC,mA))
                resize!(nzvalC, nnzC + max(nnzC,mA))
                nnzC = length(nzvalC)
            end
            colptrC[i] = ip
            for jp in colptrB[i]:(colptrB[i+1] - 1)
                nzB = nzvalB[jp]
                j = rowvalB[jp]
                for kp in colptrA[j]:(colptrA[j+1] - 1)
                    nzC = nzvalA[kp] * nzB
                    k = rowvalA[kp]
                    if xb[k] != i
                        rowvalC[ip] = k
                        ip += 1
                        xb[k] = i
                        x[k] = nzC
                    else
                        x[k] += nzC
                    end
                end
            end
            idx = colptrC[i]
            for vp in colptrC[i]:(ip - 1)
                col = rowvalC[vp]
                if col == i || x[col] > eps
                  nzvalC[idx] = x[col]
                  rowvalC[idx] = col
                  idx += 1
                end
            end
            ip = idx
        end
        colptrC[nB+1] = ip
    end

    deleteat!(rowvalC, colptrC[end]:length(rowvalC))
    deleteat!(nzvalC, colptrC[end]:length(nzvalC))

    # The Gustavson algorithm does not guarantee the product to have sorted row indices.
    Cunsorted = SparseMatrixCSC(mA, nB, colptrC, rowvalC, nzvalC)
    C = Base.SparseMatrix.sortSparseMatrixCSC!(Cunsorted, sortindices=sortindices)
    return C
end

function gershgorin(A :: SparseMatrixCSC)
  hsize = size(A, 1)
  eMin = 10000
  eMax = -10000

  for i=1:hsize
    sumM = 0
    sumP = 0
    for j=A.colptr[i]:A.colptr[i+1]-1
      hx = abs(A.nzval[j])
      sumM += hx
      if A.rowval[j] == i
        sumP = A.nzval[j]
        sumM -= hx
      end
    end
    eMax = max(eMax, sumP + sumM)
    eMin = min(eMin, sumP - sumM)
  end

  return eMax, eMin
end

function CoSP2_ref(X)
  m = size(X, 1)
  occ = m/2
  idemTol = 1e-14
  eps = 1e-5

  eMax, eMin = gershgorin(X')
  X = SparseMatrixCSC{Float64, Int32}((eMax*speye(m) - X)/(eMax - eMin))
  assert(issym(X))
  idempErr = 0
  idempErr1 = 0
  idempErr2 = 0
  BreakLoop = 0
  iter = 0
  trXOLD = 0

  t0 = time()
  spgemm_time = 0
  spadd_time = 0

  while BreakLoop == 0
    trX = trace(X)

    spgemm_time -= time()
    #X2 = X*X'
    X2 = spmatmul_witheps(X, X', eps) # CoSP2 does approximate spgemm
    spgemm_time += time()

    trX2 = trace(X2)
    trXOLD = trX
    tempX = X

    limDiff = abs(trX2 - occ) - abs(2*trX - trX2 - occ)
    if limDiff > idemTol

      spadd_time -= time()
      X = 2*X - X2
      spadd_time += time()

      trX = 2*trX - trX2
    elseif limDiff < -idemTol
      X = X2 # output X is dead here, so we can free
      trX = trX2
    else
      trX = trXOLD
      BreakLoop = 1
    end

    idempErr2 = idempErr1
    idempErr1 = idempErr
    idempErr = abs(trX - trXOLD)
    #println(idempErr)

    iter += 1

    if iter >= 25 && idempErr >= idempErr2
      BreakLoop = 1
    end
  end

  X = 2*X

  println("D Sparsity AAN = $(sum(X)), fraction = $(sum(X)/(m*m)) avg = $(sum(X)/m), max = $(maximum(X))")
  println("Number of iterations = $iter")
  println("Total time $(time() - t0), SpGEMM time $(spgemm_time), SpADD time $(spadd_time)")
end

function CoSP2_call_replacement(X)
  m = size(X, 1)
  occ = m/2
  idemTol = 1e-14
  eps = 1e-5

  eMax, eMin = gershgorin(X')
  X = SparseMatrixCSC{Float64, Int32}((eMax*speye(m) - X)/(eMax - eMin))
  assert(issym(X))
  idempErr = 0
  idempErr1 = 0
  idempErr2 = 0
  BreakLoop = 0
  iter = 0
  trXOLD = 0

  t0 = time()
  spgemm_time = 0
  spadd_time = 0

  while BreakLoop == 0
    #trX = trace(X)
    trX = SparseAccelerator.trace(X)

    spgemm_time -= time()
    #X2 = X*X'
    #X2 = spmatmul_witheps(X, X', eps) # CoSP2 does approximate spgemm
    X2 = SparseAccelerator.SpSquareWithEps(X, eps)
    spgemm_time += time()

    # trX2 = trace(X2)
    trX2 = SparseAccelerator.trace(X2)
    trXOLD = trX
    tempX = X

    limDiff = abs(trX2 - occ) - abs(2*trX - trX2 - occ)
    if limDiff > idemTol

      spadd_time -= time()
      #X = 2*X - X2 # input X is dead here, so we can free
      X = SparseAccelerator.SpAdd(2, X, -1, X2)
      spadd_time += time()

      trX = 2*trX - trX2
    elseif limDiff < -idemTol
      X = X2 # output X is dead here, so we can free
      trX = trX2
    else
      trX = trXOLD
      BreakLoop = 1
    end

    idempErr2 = idempErr1
    idempErr1 = idempErr
    idempErr = abs(trX - trXOLD)
    #println(idempErr)

    iter += 1

    if iter >= 25 && idempErr >= idempErr2
      BreakLoop = 1
    end
  end

  X = 2*X

  println("D Sparsity AAN = $(sum(X)), fraction = $(sum(X)/(m*m)) avg = $(sum(X)/m), max = $(maximum(X))")
  println("Number of iterations = $iter")
  println("Total time $(time() - t0), SpGEMM time $(spgemm_time), SpADD time $(spadd_time)")
end

function CoSP2_call_replacement_and_context_opt(X)
  m = size(X, 1)
  occ = m/2
  idemTol = 1e-14
  eps = 1e-5

  eMax, eMin = gershgorin(X')
  X = SparseMatrixCSC{Float64, Int32}((eMax*speye(m) - X)/(eMax - eMin))
  assert(issym(X))
  idempErr = 0
  idempErr1 = 0
  idempErr2 = 0
  BreakLoop = 0
  iter = 0
  trXOLD = 0

  t0 = time()
  spgemm_time = 0
  spadd_time = 0

  mknobX = (SparseAccelerator.new_matrix_knob)(X, false, false, true, true, false, false) # X is symmetric

  fknob_spgemm = (SparseAccelerator.new_function_knob)()
  (SparseAccelerator.add_mknob_to_fknob)(mknobX, fknob_spgemm)

  while BreakLoop == 0
    #trX = trace(X)
    trX = SparseAccelerator.trace(X)

    spgemm_time -= time()
    #X2 = X*X'
    #X2 = spmatmul_witheps(X, X', eps) # CoSP2 does approximate spgemm
    X2 = SparseAccelerator.SpSquareWithEps(X, eps, fknob_spgemm)
    spgemm_time += time()

    # trX2 = trace(X2)
    trX2 = SparseAccelerator.trace(X2)
    trXOLD = trX
    tempX = X

    limDiff = abs(trX2 - occ) - abs(2*trX - trX2 - occ)
    if limDiff > idemTol

      spadd_time -= time()
      #X = 2*X - X2 # input X is dead here, so we can free
      X = SparseAccelerator.SpAdd(2, X, -1, X2)
      spadd_time += time()

      trX = 2*trX - trX2
    elseif limDiff < -idemTol
      X = X2 # output X is dead here, so we can free
      trX = trX2
    else
      trX = trXOLD
      BreakLoop = 1
    end

    idempErr2 = idempErr1
    idempErr1 = idempErr
    idempErr = abs(trX - trXOLD)
    #println(idempErr)

    iter += 1

    if iter >= 25 && idempErr >= idempErr2
      BreakLoop = 1
    end
  end

  X = 2*X

  println("D Sparsity AAN = $(sum(X)), fraction = $(sum(X)/(m*m)) avg = $(sum(X)/m), max = $(maximum(X))")
  println("Number of iterations = $iter")
  println("Total time $(time() - t0), SpGEMM time $(spgemm_time), SpADD time $(spadd_time)")
end

file_name = "hmatrix.1024.mtx"
if length(ARGS) >= 1
  file_name = ARGS[1]
end

X = matrix_market_read(file_name)
assert(issym(X))

# currently, something wrong with spmatmul_witheps used in CoSP2_ref.
# So, ignore the results of CoSP2_ref
println("\nOriginal:")
CoSP2_ref(X)
CoSP2_ref(X)
println("End original.")

# Expected results: D Sparsity AAN = 12212.785128790038, fraction = 8.088207992441149e-5 avg = 0.9938789981111684, max = 1.2808837088549991
# Number of iterations = 25
println("\nCoSP2_call_replacement:")
CoSP2_call_replacement(X)
CoSP2_call_replacement(X)
println("End CoSP2_call_replacement.")

println("\nCoSP2_call_replacement_and_context_opt:")
CoSP2_call_replacement_and_context_opt(X)
CoSP2_call_replacement_and_context_opt(X)
println("End CoSP2_call_replacement_and_context_opt.")
