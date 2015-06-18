using DataFrames
using Compose
using Gadfly

const LIB_PATH = "/media/sf_VBoxVMShared/SparseAccelerator/lib/libcsr.so"

function CSR_ReorderMatrix(A::SparseMatrixCSC, newA::SparseMatrixCSC, P::Vector, Pprime::Vector, getPermutation::Bool, oneBasedInput::Bool, oneBasedOutput::Bool)
  ccall((:CSR_ReorderMatrix, LIB_PATH), Void,
              (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
               Ptr{Cint}, Ptr{Cint}, Bool, Bool, Bool),
               A.m, A.n, pointer(A.colptr), pointer(A.rowval), pointer(A.nzval),
               pointer(newA.colptr), pointer(newA.rowval), pointer(newA.nzval),
               pointer(P), pointer(Pprime), getPermutation, oneBasedInput, oneBasedOutput)
end


N = 14
DENSITY = 0.2 #0.01

A = sprand(N, N, DENSITY)

# align with library's assumptions: type and symmetry
A=convert(SparseMatrixCSC{Cdouble, Cint}, A)
for i=1:N for j=1:N A[i, j] = A[j, i] end end

P = Array(Cint,size(A,2))
Pprime = Array(Cint,size(A,2))
Aprime = SparseMatrixCSC(A.m,A.n,Array(Cint,size(A.colptr,1)),Array(Cint,size(A.rowval,1)),Array(Cdouble,size(A.nzval,1)))
CSR_ReorderMatrix(A,Aprime,P,Pprime,true,true,true)

#println("P = ", P)
#println("A = ", A)
#println("Aprime = ", Aprime)

is, js, vs = findnz(A)
d = DataFrame(i=is, j=js, v=vs)
p=plot(d, x="j", y="i", color="v",
			     Coord.cartesian(yflip=true),
			     Scale.color_continuous,
			     Scale.x_continuous,
			     Scale.y_continuous,
			     Geom.rectbin,
			     Stat.identity)
c = render(p)	
draw(SVG("matrix0.svg", 7inch, 7inch), c)
run(`convert matrix0.svg matrix0.png`)
	    
is, js, vs = findnz(Aprime)
d = DataFrame(i=is, j=js, v=vs)
p=plot(d, x="j", y="i", color="v",
			 Coord.cartesian(yflip=true),
			 Scale.color_continuous,
			 Scale.x_continuous,
			 Scale.y_continuous,
			 Geom.rectbin,
			 Stat.identity)
c = render(p)	
draw(SVG("Aprime.svg", 7inch, 7inch), c)
run(`convert Aprime.svg Aprime.png`)

processed = Vector{Bool}(N)
for i = 1:N processed[i] = false end

# Let us say the permutation (P) is:
#    1 2 3 4 5
#    5 3 4 1 2
# Then we expect that Row(column) 4 appears in row(column) 1, 5 appears in 2, 2 in 3, 3 in 4, and 1 in 5.
# To do that, we should visit the nodes in the order of 1->4->3->2->5. That is from Pprime.
   
IDX = 1
while true
	found = false
	start = 1
	for i = 1:N
		if !processed[i]
			found = true
			start = i
			break
		end
	end
	if !found
		break
	end
	
	i = start
	while Pprime[i] + 1 != start
	    # generate a permuation matrix
	    E = eye(N, N)
	    E[i, i] = 0
	    E[i, Pprime[i] + 1] = 1

	    E[Pprime[i] + 1, Pprime[i] + 1] = 0			
	    E[Pprime[i] + 1, i] = 1

	    # permutate rows and columns
	    A = E * A * E
#	    println("A after itr ", IDX, "(", i, "<->", Pprime[i]+1, "): ", A)

	    is, js, vs = findnz(A)
	    d = DataFrame(i=is, j=js, v=vs)
	    p=plot(d, x="j", y="i", color="v",
			     Coord.cartesian(yflip=true),
			     Scale.color_continuous,
			     Scale.x_continuous,
			     Scale.y_continuous,
			     Geom.rectbin,
			     Stat.identity)
	    c = render(p)	
	    draw(SVG("matrix$IDX.svg", 7inch, 7inch), c)
	    run(`convert matrix$IDX.svg matrix$IDX.png`)

	    processed[i] = true
	    IDX +=1	  
	    i = Pprime[i] + 1
	end
	processed[i] = true	    
end	

#convert -delay 0.01 -loop 0 'matrix%d.png[0-100]' animation.gif
run(`convert -delay 0.01 -loop 0 \'matrix%d.png[0:N]\' animation.gif`)
