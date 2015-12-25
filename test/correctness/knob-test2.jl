# julia knob-test2.jl ../matrices/bcsstk14.mtx
# Turning on any of the commented statement to repro the bug.
 
include("../../src/SparseAccelerator.jl")
include("./utils.jl")
using SparseAccelerator

function pcg_symgs(x, A, b, tolerance, maxiter)
    __mknobL4__ = SparseAccelerator.new_matrix_knob(:L,true,true,false,false,false,true)
    __mknobA1__ = SparseAccelerator.new_matrix_knob(:A,true,true,true,true,false,true)
    __mknobU6__ = SparseAccelerator.new_matrix_knob(:U,true,true,false,false,false,true)
    __fknobilu__ = SparseAccelerator.new_function_knob()
    __fknob3__ = SparseAccelerator.new_function_knob()
    __fknob2__ = SparseAccelerator.new_function_knob()
    __fknob5__ = SparseAccelerator.new_function_knob()
    __fknob7__ = SparseAccelerator.new_function_knob()
    __fknob0__ = SparseAccelerator.new_function_knob()
    __fknob8__ = SparseAccelerator.new_function_knob()
    SparseAccelerator.set_derivative(__mknobL4__,SparseAccelerator.DERIVATIVE_TYPE_LOWER_TRIANGULAR,__mknobA1__)
    SparseAccelerator.set_derivative(__mknobU6__,SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC,__mknobA1__)
    SparseAccelerator.set_derivative(__mknobL4__,SparseAccelerator.DERIVATIVE_TYPE_SYMMETRIC,__mknobA1__)
    SparseAccelerator.set_derivative(__mknobU6__,SparseAccelerator.DERIVATIVE_TYPE_UPPER_TRIANGULAR,__mknobA1__)
    SparseAccelerator.add_mknob_to_fknob(__mknobL4__,__fknob3__)
    SparseAccelerator.add_mknob_to_fknob(__mknobA1__,__fknob2__)
    SparseAccelerator.add_mknob_to_fknob(__mknobU6__,__fknob5__)
    SparseAccelerator.add_mknob_to_fknob(__mknobL4__,__fknob7__)
    SparseAccelerator.add_mknob_to_fknob(__mknobA1__,__fknob0__)
    SparseAccelerator.add_mknob_to_fknob(__mknobU6__,__fknob8__)
    SparseAccelerator.add_mknob_to_fknob(__mknobA1__,__fknobilu__)
    total_time = time()
    
    
    L = tril(A)
    U = SparseMatrixCSC{Cdouble, Cint}(spdiagm(1./diag(A)))*triu(A)
    # ISSUE: Either of the following will make the code run endlessly
    # L, U = SparseAccelerator.ilu(A)
    # L, U = SparseAccelerator.ilu(A, __fknobilu__)
    
    r = b - SparseAccelerator.SpMV(1,A,x,__fknob2__)
    normr0 = SparseAccelerator.norm(r)
    z = copy(r)
    SparseAccelerator.fwdTriSolve!(z,L,r,__fknob3__)
    SparseAccelerator.bwdTriSolve!(z,U,z,__fknob5__)
    p = copy(z)
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
            SparseAccelerator.fwdTriSolve!(z,L,r,__fknob7__)
            SparseAccelerator.bwdTriSolve!(z,U,z,__fknob8__)
            rz = SparseAccelerator.dot(r,z)
            SparseAccelerator.WAXPBY!(p,1,z,rz / old_rz,p)
            k = k + 1
    end
    total_time = time() - total_time
    println("total = ",total_time,"s")
    SparseAccelerator.delete_matrix_knob(__mknobL4__)
    SparseAccelerator.delete_matrix_knob(__mknobA1__)
    SparseAccelerator.delete_matrix_knob(__mknobU6__)
    SparseAccelerator.delete_function_knob(__fknob3__)
    SparseAccelerator.delete_function_knob(__fknob2__)
    SparseAccelerator.delete_function_knob(__fknob5__)
    SparseAccelerator.delete_function_knob(__fknob7__)
    SparseAccelerator.delete_function_knob(__fknob0__)
    SparseAccelerator.delete_function_knob(__fknob8__)
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
println("\tOptimized sum of x=", sum(x))
println("\tOptimized k=", k)
println("\tOptimized rel_err=", rel_err)