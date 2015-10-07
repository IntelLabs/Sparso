#!/bin/bash

. /opt/intel/compiler/latest/bin/compilervars.sh # make sure you also did this before compiling libcsr.so 

# Use 1 socket, 1 HW thread per core
export OMP_NUM_THREADS=14
export KMP_AFFINITY=granularity=fine,compact,1

for i in `cat ipm.lst`; do
  echo !!!!$i
  #echo %%%%manual
  #julia context-test4.jl ~/matrices/lp/$i/$i # run this if you want to compare with manually optimized version
  #echo
  #echo %%%%auto
  julia context-test3.jl ~/matrices/lp/$i/$i
  echo
  echo
done
