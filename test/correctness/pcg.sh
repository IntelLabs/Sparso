#!/bin/bash

. /opt/intel/compiler/latest/bin/compilervars.sh intel64 # make sure you also did this before compiling libcsr.so 
# Use 1 socket, 1 HW thread per core
export OMP_NUM_THREADS=14
export KMP_AFFINITY=granularity=fine,compact,1

for i in `cat pcg.lst`; do
  echo !!!!$i
  echo %%%%auto-without-context-serial
  OMP_NUM_THREADS=1 julia context-test5-without-context.jl ~/matrices/$i.mtx
  echo %%%%auto-without-context
  julia context-test5-without-context.jl ~/matrices/$i.mtx
  echo
  echo %%%%auto-without-reordering
  julia context-test5-without-reordering.jl ~/matrices/$i.mtx
  echo
  echo %%%%auto
  julia context-test5.jl ~/matrices/$i.mtx
  echo
  echo
done
