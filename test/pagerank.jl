using MatrixMarket

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

A = MatrixMarket.mmread(ASCIIString(ARGS[1]))
A = spones(A)

m = size(A, 1)
p = repmat([1/m], m)
@time pagerank(A, p, 0.15)
