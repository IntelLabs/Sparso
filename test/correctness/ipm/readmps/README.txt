This directory contains MATLAB code for reading files in MPS format.  The
code supports integer and linear programming problems as well as quadratic
programming problems.  
 
In order to use this routine you must understand the very general form of
a mixed integer linear/quadratic program used in the MPS format.  This is 
 
min/max (1/2) x'*Q*x + obj'*x 
                 b-r (<=, =, >=) Ax  (<=, =, >=) b+r
                 l<=x<=u

Rows and columns of A are given names.  A major function of this code is to
maintain tables that associate these names with row and column numbers.  
Rows can represent equality constraints, less than or equal to constraints,
greater than or equal to constraints, or (in conjunctions with ranges) the
row value can be constrained to lie within a range.  

There may be multiple objective functions, multiple right hand sides, 
multiple range vectors, and multiple upper and lower bound vectors in
an MPS file.  The user should have some way outside of the MPS file to 
select which of the many objective functions, right hand side, ranges, 
and bounds are to be used in solving the problem.   

The useage is simply 
 
problem=readmps(filename)

or (for files in fixed format that might have spaces in row/column names)

problem=readmpsfx(filename)

Fields of the problem output
 
name                            Problem name.
objsesense                      'MINIMIZE', 'MIN', 'MAX', or 'MAXIMIZE'. 
objname                         Name of the objective function row.
problem.refrow                  Name of the reference row for SOS's.
rownames                        Cell array of row names
rowtypes                        Cell array of row types ('L','G','N','E')
columnnames                     Cell array of column names
boundnames                      Cell array of names of bounds
rhsnames                        Cell array of names of right hand sides
rangenames                      Cell array of names of ranges
lbnds                           Sparse array of lower bounds
ubnds                           Sparse array of upper bounds
rhs                             Sparse array of right hand sides
ranges                          Sparse array of ranges
bintflags                       Sparse array of flags.  bintflags(i,j)=1 if
                                column j in bound set i is an integer column
                                in bound set i (different bound sets might have
                                different integer columns.)
intflags                        intflags(j)=1 if column j is an integer column
sos1flags                       sos1flags(j)=1 if column j is in an SOS1
sos2flags                       sos1flags(j)=1 if column j is in an SOS2
sos3flags                       sos1flags(j)=1 if column j is in an SOS3
Q                               Sparse array of quadratic objective function
                                coefficients
rowtable                        hash table for row names
coltable                        hash table for column names
boundtable                      hash table for bound names
rhstable                        hash table for right hand side names
rangetable                      hash table for range names

The hash tables can be used with the tablelookup() function to look up
a row (column, bound set, right hand side, or range) name and return 
the corresponding row (etc.) number.

Update History:

Version 2, 1/18/2015.  This version uses the first "N" row as the default 
objective row if none is specified.  The hash function used with the hash
tables has been improved to reduce frequent collisions that were happening on 
certain problems.  Added support for fixed format files.  
