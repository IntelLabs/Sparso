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

include("../../../src/Sparso.jl")
using Sparso

function spmatmul_witheps_and_flops{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}, eps;
                         sortindices::Symbol = :sortcols)
    mA, nA = size(A)
    mB, nB = size(B)
    nA==mB || throw(DimensionMismatch())

    flops = 0

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

                    flops += 2
                end
            end
            idx = colptrC[i]
            for vp in colptrC[i]:(ip - 1)
                col = rowvalC[vp]
                if col == i || abs(x[col]) > eps
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
    return C, flops
end

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
                if col == i || abs(x[col]) > eps
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

@doc """ Using Gershgorin disc to estimate the bounds of eigen value. """
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

function count_CoSP2_flop(X)
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

  flops = 0
  spgemm_flops = 0

  while BreakLoop == 0
    trX = trace(X)
    flops += m

    X2, f = spmatmul_witheps_and_flops(X, X', eps) # CoSP2 does approximate spgemm
    flops += f
    spgemm_flops += f

    trX2 = trace(X2)
    flops += m
    trXOLD = trX

    limDiff = abs(trX2 - occ) - abs(2*trX - trX2 - occ)
    if limDiff > idemTol

      X = 2*X - X2
      flops += nnz(X) + nnz(X2)

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
  flops += nnz(X)

  println("X sum = $(sum(X)), max = $(maximum(X))")
  println("Number of iterations = $iter")
  #println("$flops $spgemm_flops")
  return flops
end

function CoSP2_ref(X, flops)
  set_matrix_property(Dict(
      :X => SA_SYMM_VALUED, 
    )
  )

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

  total_time = -time()
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

  total_time += time()

  println("X sum = $(sum(X)), max = $(maximum(X))")
  println("Number of iterations = $iter")
  println("Total time $total_time ($(flops/total_time/1e9) gflops), SpGEMM time $spgemm_time, SpADD time $spadd_time")
end

file_name = "hmatrix.512.mtx"
if length(ARGS) >= 1
  file_name = ARGS[1]
end

X = matrix_market_read(file_name)
assert(issym(X))

flops = count_CoSP2_flop(X)

if length(ARGS) == 2
  test = ARGS[2]
else
  test = "julia"
end

if test == "call-repl"
  set_options(SA_ENABLE, SA_USE_SPMP, SA_REPLACE_CALLS)
elseif test == "context"
  set_options(SA_ENABLE, SA_USE_SPMP, SA_CONTEXT, SA_REPLACE_CALLS)
elseif test == "reorder"
  set_options(SA_ENABLE, SA_USE_SPMP, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)
elseif test == "verbose"
  set_options(SA_ENABLE, SA_USE_SPMP, SA_VERBOSE, SA_CONTEXT, SA_REORDER, SA_REPLACE_CALLS)
end

println("compiler warmup (ignored): ")
if test == "julia"
  CoSP2_ref(X, flops)
else
  @acc CoSP2_ref(X, flops)
end

println("\nRUN: ")
if test == "julia"
  CoSP2_ref(X, flops)
else
  @acc CoSP2_ref(X, flops)
end
