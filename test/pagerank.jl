using MatrixMarket

function mmread(filename::ASCIIString, infoonly::Bool)
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
    mmfile = open(filename,"r")
    tokens = split(chomp(readline(mmfile)))
    if length(tokens) != 5 error("Not enough words on header line") end
    if tokens[1] != "%%MatrixMarket" error("Not a valid MatrixMarket header.") end
    (head1, rep, field, symm) = map(lowercase, tokens[2:5])
    if head1 != "matrix"
        error("This seems to be a MatrixMarket $head1 file, not a MatrixMarket matrix file")
    end
    if field != "real" && field != "pattern" error("complex fields not yet allowed") end

    ll   = readline(mmfile)         # Read through comments, ignoring them
    while length(ll) > 0 && ll[1] == '%' ll = readline(mmfile) end
    dd     = int(split(ll))         # Read dimensions
    rows   = dd[1]
    cols   = dd[2]
    entries = rep == "coordinate" ? dd[3] : rows * cols
    if infoonly return rows, cols, entries, rep, field, symm end
    if rep == "coordinate"
        rr = Array(Int, entries)
        cc = Array(Int, entries)
        xx = Array(Float64, entries)
        for i in 1:entries
            flds = split(readline(mmfile))
            rr[i] = int32(flds[1])
            cc[i] = int32(flds[2])
            if field != "pattern"
              xx[i] = float64(flds[3])
            else
              xx[i] = 1
            end
        end
        return sparse(rr, cc, xx, rows, cols)
    end
    reshape([float64(readline(mmfile)) for i in 1:entries], (rows,cols))
end

mmread(filename::ASCIIString) = mmread(filename, false)

function pagerank(A, p, r) # p: initial rank, r: damping factor
  t = 0
  t2 = 0
  t3 = 0
  t4 = 0

  d = vec(sum(A, 2)) # num of neighbors

  for i = 1:100
    q = p./d

    tic()
    Aq = A*q
    t += toq()

    p2 = r + (1-r)*Aq
    println("$i: $(norm(p - p2)/norm(p))") # print out convergence

    p = p2
  end

  println("SpMV takes $t sec.")
end

A = mmread(ASCIIString(ARGS[1]))
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
@time pagerank(A, p, 0.15)
