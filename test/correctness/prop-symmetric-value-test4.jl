include("../../src/SparseAccelerator.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT)

include("utils.jl")

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

function foo()
    eps = 1e-5

    m = 10
    X = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    Y = generate_symmetric_nonzerodiagonal_sparse_matrix(m)
    # X = SparseMatrixCSC{Float64, Int32}((eMax*speye(m) - X)/(eMax - eMin))

    for i = 1:2
        A = spmatmul_witheps(X, X', eps)
        B = spmatmul_witheps(X', X, eps)
        C = spmatmul_witheps(X, Y', eps)
        D = spmatmul_witheps(X, Y, eps)
    end
end

@acc foo()
