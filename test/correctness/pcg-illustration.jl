# Note: tested with julia-b77587362d (0.5.0 dev) nightly build. Julia-0.4 has
# issues with "\", and cannot run it through.
# Run with
#   julia-b77587362d/bin/julia pcg-illustration.jl  ../matrices/bcsstk14.mtx
# Look at the comments for issues.

include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

set_options(SA_ENABLE, SA_VERBOSE, SA_USE_SPMP, SA_CONTEXT, SA_CONTEXT_FUNC, SA_REPLACE_CALLS, SA_REORDER)

function pcg_symgs(x, A, b, tolerance, maxiter)
    set_matrix_property(:L, SA_LOWER_OF, :A)
    set_matrix_property(:U, SA_UPPER_OF, :A)
    set_matrix_property(Dict(
        :A => SA_CONST_VALUED | SA_SYMM_STRUCTURED | SA_SYMM_VALUED,
        :L => SA_CONST_VALUED,
        :U => SA_CONST_VALUED,
        )
    )

    total_time = time()

    # TODO: replace with ILU0 or ICHOL
    L = tril(A)
    U = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)
    r = b - A * x
    normr0 = norm(r)
    z = copy(r) # A workaround when SA_CONTEXT_FUNC is enabled so that z is defined before the following statement is replaced with fwdTriSolve!(z, ..)
    z = L \ r    
    z = U \ z    
    p = copy(z)
    rz = dot(r, z)

    rel_err = 1
    k = 1
    while true
        old_rz = rz
        Ap = A*p
        alpha = old_rz / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rel_err = norm(r)/normr0
        if rel_err < tolerance 
            break
        end
        z = L \ r
        z = U \ z
        rz = dot(r, z)
        p = z + (rz/old_rz) * p
        k += 1
    end

    total_time = time() - total_time
    println("total = $(total_time)s")

    return x, k, rel_err
end

function pcg_symgs_opt(x, A, b, tolerance, maxiter)
    __mknobL3__ = SparseAccelerator.new_matrix_knob(:L,true,true,false,false,false,false)
    __mknobU5__ = SparseAccelerator.new_matrix_knob(:U,true,true,false,false,false,false)
    __mknobA1__ = SparseAccelerator.new_matrix_knob(:A,true,true,true,true,false,false)
    __fknob4__ = SparseAccelerator.new_function_knob()
    __fknob2__ = SparseAccelerator.new_function_knob()
    __fknob0__ = SparseAccelerator.new_function_knob()
    __fknob33__ = SparseAccelerator.new_function_knob()
    __fknob22__ = SparseAccelerator.new_function_knob()
    __fknob55__ = SparseAccelerator.new_function_knob()
    SparseAccelerator.add_mknob_to_fknob(__mknobU5__,__fknob4__)
    SparseAccelerator.add_mknob_to_fknob(__mknobL3__,__fknob2__)
    SparseAccelerator.add_mknob_to_fknob(__mknobA1__,__fknob0__)
    SparseAccelerator.add_mknob_to_fknob(__mknobL3__,__fknob33__)
    SparseAccelerator.add_mknob_to_fknob(__mknobA1__,__fknob22__)
    SparseAccelerator.add_mknob_to_fknob(__mknobU5__,__fknob55__)
    SparseAccelerator.set_derivative(__mknobL3__,SparseAccelerator.DERIVATIVE_TYPE_LOWER_TRIANGULAR,__mknobA1__)
    SparseAccelerator.set_derivative(__mknobU5__,SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC,__mknobA1__)
    SparseAccelerator.set_derivative(__mknobL3__,SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC,__mknobA1__)
    SparseAccelerator.set_derivative(__mknobU5__,SparseAccelerator.DERIVATIVE_TYPE_UPPER_TRIANGULAR,__mknobA1__)
    SparseAccelerator.set_reordering_decision_maker(__fknob2__, SparseAccelerator.PERM_VECTORS_ARE_INVERSE)
    __reordering_status6__ = [false,Main.C_NULL,Main.C_NULL,Main.C_NULL,Main.C_NULL,0.0]
    total_time = Main.time()
    L = Main.tril(A)
    U = SparseMatrixCSC{Cdouble,Cint}(Main.spdiagm(1 ./ Main.diag(A))) * Main.triu(A)
    r = b - SparseAccelerator.SpMV(1,A,x,__fknob22__)
    normr0 = SparseAccelerator.norm(r)
    z = Main.copy(r)
    
    # ISSUE: turning on the following two, and comment out the next two, statements,
    # will cause the code to run endlessly. It may be because fwdTriSolve! has cached
    # something at this time. However, later, in the loop, fwdTriSolve! may cache
    # other things during reordering, and the informations are not consistent with
    # each other. This is my guess.
#    SparseAccelerator.fwdTriSolve!(z,L,r,__fknob33__)
#    SparseAccelerator.bwdTriSolve!(z,U,z,__fknob55__)
    z = L \ r
    z = U \ z

    p = Main.copy(z)
    rz = SparseAccelerator.dot(r,z)
    rel_err = 1
    k = 1
    while true
            old_rz = rz
            Ap = SparseAccelerator.SpMV(1,A,p,__fknob0__)
            alpha = old_rz / SparseAccelerator.dot(p,Ap)
            SparseAccelerator.WAXPBY!(x,1,x,alpha,p)
            SparseAccelerator.WAXPBY!(r,1,r,-alpha,Ap)
            rel_err = SparseAccelerator.norm(r) / normr0
            if rel_err < tolerance
                break
            end
            SparseAccelerator.fwdTriSolve!(z,L,r,__fknob2__)
            SparseAccelerator.reordering(__fknob2__,__reordering_status6__,A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobA1__,U,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobU5__,:__delimitor__,p,SparseAccelerator.ROW_PERM,r,SparseAccelerator.ROW_PERM,x,SparseAccelerator.ROW_PERM)
            
            # ISSUE: turning on this statement will cause the code to run endlessly. The only difference from the above statement is in U's permutation vectors. 
            # However, this statement is actually right, because from the source code: "z=L\r" and "dot(r,z)", we know "Lz=r" and "dot(r,z)",
            # and therefore: 
            #    L' row permutation vector is the same as r's, which is the same as z's, which is inverse to L's column permutation vector.
            # Therefore,
            #    L's row and column permutation vector must be inverse to each other. 
            # However, this is the output of when LOG_REORDERING in knob.cpp is true:
            #    TriangularSolve: row_perm=0x3966d40 row_inverse_perm=0x39c7e00 col_perm=0x3966d40 col_inverse_perm=0x39c7e00
            # Obviously, the row and column perm vectors are the same, instead of inverse to each other. 
            # TODO: I added an option perm_restriction to set_reordering_decision_maker(). The option has the combination of 3 possible
            # values: "PERM_VECTORS_ARE_INVERSE" or "PERM_VECTORS_ARE_SAME" or "PERM_VECTORS_ARE_INDEPENDENT". Please see how
            # you can utilize this info.
            # Note: I guess usually, just one value is possible, but not sure if more than one restriction exists, for example:
            #    "PERM_VECTORS_ARE_INVERSE" and "PERM_VECTORS_ARE_SAME"
            # There are permutations that can satisfy both restrictions, like identity permutation matrix.
#            SparseAccelerator.reordering(__fknob2__,__reordering_status6__,A,SparseAccelerator.ROW_PERM,SparseAccelerator.ROW_INV_PERM,__mknobA1__,U,SparseAccelerator.COL_INV_PERM,SparseAccelerator.COL_PERM,__mknobU5__,:__delimitor__,p,SparseAccelerator.ROW_PERM,r,SparseAccelerator.ROW_PERM,x,SparseAccelerator.ROW_PERM)

            SparseAccelerator.bwdTriSolve!(z,U,z,__fknob4__)
            rz = SparseAccelerator.dot(r,z)
            SparseAccelerator.WAXPBY!(p,1,z,rz / old_rz,p)
            k = k + 1
    end
    total_time = Main.time() - total_time
    Main.println("total = ",total_time,"s")
    SparseAccelerator.delete_matrix_knob(__mknobL3__)
    SparseAccelerator.delete_matrix_knob(__mknobU5__)
    SparseAccelerator.delete_matrix_knob(__mknobA1__)
    SparseAccelerator.delete_function_knob(__fknob4__)
    SparseAccelerator.delete_function_knob(__fknob2__)
    SparseAccelerator.delete_function_knob(__fknob0__)
    SparseAccelerator.reverse_reordering(__reordering_status6__,:__delimitor__,x,SparseAccelerator.ROW_PERM)
    return x,k,rel_err
end

A         = matrix_market_read(ARGS[1], true, true)
m         = size(A, 1)
b         = ones(Float64, m)
originalA = copy(A)
tolerance = 1e-7
maxiter   = 20000

x = zeros(Float64, m)
x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
println("\tOriginal sum of x=", sum(x))
println("\tOriginal k=", k)
println("\tOriginal rel_err=", rel_err)

x = zeros(Float64, m)
b = ones(Float64, m)
A = copy(originalA)
x, k, rel_err = pcg_symgs_opt(x, A, b, tolerance, maxiter)
println("\tManual opt version: sum of x=", sum(x))
println("\tManual opt version:  k=", k)
println("\tManual opt version: rel_err=", rel_err)
flush(STDOUT)

# This will run endlessly. Please look at the comments in pcg_symgs_opt() and
# play with it. Once the issues are solved there, this should work automatically.
#x = zeros(Float64, m)
#b = ones(Float64, m)
#A = copy(originalA)
#@acc x, k, rel_err = pcg_symgs(x, A, b, tolerance, maxiter)
#println("\tOptimized version: sum of x=", sum(x))
#println("\tOptimized version:  k=", k)
#println("\tOptimized version: rel_err=", rel_err)
