module MatrixMarket

function mmread(filename, infoonly::Bool)
#      Reads the contents of the Matrix Market file 'filename'
#      into a matrix, which will be either sparse or dense,
#      depending on the Matrix Market format indicated by
#      'coordinate' (coordinate sparse storage), or
#      'array' (dense array storage).  The data will be duplicated
#      as appropriate if symmetry is indicated in the header. (Not yet
#      implemented).
#
#      If infoonly is true information on the size and structure is
#      returned.
    if infoonly
        filename1 = ASCIIString(filename)
        mmfile = open(filename1,"r")
        tokens = split(chomp(readline(mmfile)))
        if length(tokens) != 5 error("Not enough words on header line") end
        if tokens[1] != "%%MatrixMarket" error("Not a valid MatrixMarket header.") end
        (head1, rep, field, symm) = map(lowercase, tokens[2:5])
        if head1 != "matrix"
            error("This seems to be a MatrixMarket $head1 file, not a MatrixMarket matrix file")
        end
        if field != "real" error("non-float fields not yet allowed") end
    
        ll   = readline(mmfile)         # Read through comments, ignoring them
        while length(ll) > 0 && ll[1] == '%' ll = readline(mmfile) end
        dd     = int(split(ll))         # Read dimensions
        rows   = dd[1]
        cols   = dd[2]
        entries = rep == "coordinate" ? dd[3] : rows * cols
        
        return rows, cols, entries, rep, field, symm
    end
    
    # Read into a COO array (stored statically insded mm_read)
    sizes::Vector{Cint} = [0, 0, 0, 0]
    ccall((:load_matrix_market_step1, "../lib/libcsr.so"), Void, 
        (Ptr{Uint8}, Ptr{Cint}), filename, pointer(sizes))

    is_symmetric::Cint = sizes[1] # 0/1: true/false    
    if is_symmetric == 0
        throw("mmread: matrix market file is NOT symmetric. Sparse cannot handle it now, since we need a matrix to be symmetric, and thus we can use its CSC representation as CSR")
    end
    
    m::Cint = sizes[2]
    n::Cint = sizes[3]
    assert(m ==n)
    nnz::Cint = sizes[4]  #Note: if symmetric, nnz includes elements in both upper and lower triangle 
    
    v = Array(Cdouble, nnz)
    i = Array(Cint, nnz)
    j = Array(Cint, n + 1)
    
    A = SparseMatrixCSC{Cdouble, Cint}(m, n, j, i, v)

    # Convert the COO array to CSR array.
    ccall((:load_matrix_market_step2, "../lib/libcsr.so"), Void, 
        (Ptr{Uint8}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Bool),
        filename, pointer(v), pointer(i), pointer(j), pointer(sizes), true)
        
    distance = div(nnz, 100); # print about 100 elements to check manually
    #ccall((:CSR_PrintSomeValues, "../lib/libcsr.so"), Void, 
        #(Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Bool),
        #m, n, pointer(j), pointer(i), pointer(v), convert(Cint, distance), true)

    A
end

function mmread_reorder(filename)
    mtx_file = open(filename, "r")
    modification_time::Float64 = mtime(mtx_file)
    close(mtx_file)
    
    reordered_filename = string(filename, ".reordered")
    if isfile(reordered_filename)
        # The file format is as follows:
        #   last modification time of the original mtx file
        #   m n #non_zeros Ti_num Tv_num colptr rowval nzval
        #   PermutationVector ReversePermutationVector
        reordered_file = open(reordered_filename, "r")
        modification_time1 = read(reordered_file, Float64)
        if (modification_time == modification_time1)
            # The reordered file has the up to date info
            m = read(reordered_file, Int)
            n = read(reordered_file, Int)
            nnz = read(reordered_file, Int64)
            
            Ti_num = read(reordered_file, Int)
            Tv_num = read(reordered_file, Int)
            dict = Dict{Int, Any}(0 => Int, 1 => Int32, 2 => Int64, 3 => Float64)
            Ti = dict[Ti_num]
            Tv = dict[Tv_num]
            
            colptr = read(reordered_file, Ti, n + 1)
            println("colptr is ", colptr)
            rowval = read(reordered_file, Ti, nnz)
            println("rowval is ", rowval)
            nzval  = read(reordered_file, Tv, nnz)
            println("nzval is ", nzval)
            A = SparseMatrixCSC(m, n, colptr, rowval, nzval)
            println("A is ", A)
            
            P      = read(reordered_file, Cint, n)
            println("P is ", P)
            
            Pprime = read(reordered_file, Cint, n)
            println("Pprime is ", Pprime)
            close(reordered_file)
            
            return A, P, Pprime
        else
            # The reordered file is out of date
            close(reordered_file)
            rm(reordered_filename)
        end
    end
    
    # Read the file to build a CSC array, and reorder it
    A          = mmread(filename)
    m::Int     = A.m
    n::Int     = A.n
    nnz::Int64 = size(A.nzval, 1)
    Tv         = eltype(A.nzval)
    Ti         = eltype(A.colptr)
    colptr     = Array(Ti, n + 1)
    rowval     = Array(Ti, nnz)
    nzval      = Array(Tv, nnz)
    newA       = SparseMatrixCSC{Tv, Ti}(m, n, colptr, rowval, nzval)
    
    P          = Array(Cint, n)
    Pprime     = Array(Cint, n)
    
    ccall((:CSR_ReorderMatrix, "../lib/libcsr.so"), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool, Bool, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(newA.colptr), pointer(newA.rowval), pointer(newA.nzval),
               pointer(P), pointer(Pprime), true, true, true)

            println("colptr is ", newA.colptr)
            println("rowval is ", newA.rowval)
            println("nzval is ", newA.nzval)

    # Store the reordered A and P, PPrime into a file
    reordered_file = open(reordered_filename, "w")
    
    dict = Dict{Any, Int}(Int => 0, Int32 => 1, Int64 => 2, Float64 => 3)
    write(reordered_file, modification_time)
    write(reordered_file, m)
    write(reordered_file, n)
    write(reordered_file, nnz)
    write(reordered_file, dict[Ti])
    write(reordered_file, dict[Tv])
    write(reordered_file, newA.colptr)
    write(reordered_file, newA.rowval)
    write(reordered_file, newA.nzval)
    write(reordered_file, P)
    write(reordered_file, Pprime)
    close(reordered_file)

    return newA, P, Pprime
end

mmread(filename) = mmread(filename, false)

end # module

