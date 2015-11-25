#!/bin/bash

#. /opt/intel/compiler/latest/bin/compilervars.sh intel64 # make sure you also did this before compiling libcsr.so 
# Use 1 socket, 1 HW thread per core
export OMP_NUM_THREADS=14
export KMP_AFFINITY=granularity=fine,compact,1

for i in `cat lbfgs.lst`; do
  echo !!!!$i
  echo %%%%parallel
  julia lbfgs.jl ~/matrices/$i.mtx
  echo
  echo %%%%parallel-new
  julia lbfgs-new.jl ~/matrices/$i.mtx
  echo
  echo %%%%serial
  OMP_NUM_THREADS=1 julia lbfgs-serial.jl ~/matrices/$i.mtx
  echo
done
