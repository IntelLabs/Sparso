#!/bin/bash

. /opt/intel/compiler/latest/bin/compilervars.sh intel64 # make sure you also did this before compiling libcsr.so 
# Use 1 socket, 1 HW thread per core
export OMP_NUM_THREADS=14
export KMP_AFFINITY=granularity=fine,compact,1

declare -a args=("12 5", "14 3", "16 1")

for arg in "${args[@]}"; do
  echo !!!!qbits=$arg
  echo %%%%manual-serial
  OMP_NUM_THREADS=1 julia adiabatic.jl $arg M
  echo
  echo %%%%manual-parallel
  julia adiabatic.jl $arg M
  echo
  echo %%%%auto
  julia adiabatic.jl $arg C
  echo
done
