#!/bin/bash

#. /opt/intel/compiler/latest/bin/compilervars.sh intel64 # make sure you also did this before compiling libcsr.so

# Use 1 socket, 1 HW thread per core
export OMP_NUM_THREADS=14
export KMP_AFFINITY=granularity=fine,compact,1

declare -a benchmarks=("cosp2" "ipm" "lbfgs" "pagerank" "pcg" "adiabatic")

echo "OMP_NUM_THREADS = ${OMP_NUM_THREADS}, KMP_AFFINITY = ${KMP_AFFINITY}"

for b in "${benchmarks[@]}"; do
    echo -e "\n==== RUNNING $b ====\n"
    cd $b
    . ${b}.sh
    cd ..
done
