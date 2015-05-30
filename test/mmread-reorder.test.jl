using MatrixMarket

A, P, PPrime = MatrixMarket.mmread_reorder(ARGS[1])
println("A  is: \n", A)
println("\nP  is: \n", P)
println("\nP' is: \n",PPrime)
