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
    
    v::Vector{Cdouble} = Vector{Cdouble}(nnz)
    i::Vector{Cint} = Vector{Cint}(nnz)
    j::Vector{Cint} = Vector{Cint}(n + 1)
    
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

mmread(filename) = mmread(filename, false)

end # module

