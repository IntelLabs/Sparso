#!/bin/bash

<<"LICENSE"
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE

#. /opt/intel/compiler/latest/bin/compilervars.sh intel64 # make sure you also did this before compiling libcsr.so

# Use 1 socket, 1 HW thread per core
if [[ -z ${OMP_NUM_THREADS+x} ]]; then
    export OMP_NUM_THREADS=14
fi
if [[ -z ${KMP_AFFINITY+x} ]]; then
    export KMP_AFFINITY=granularity=fine,compact,1
fi

declare -a benchmarks=("cosp2" "ipm" "lbfgs" "pagerank" "pcg" "adiabatic")

echo "==== Config: OMP_NUM_THREADS=${OMP_NUM_THREADS}, KMP_AFFINITY=${KMP_AFFINITY}"

for b in "${benchmarks[@]}"; do
    echo -e "\n==== RUNNING $b ====\n"
    cd $b
    . ${b}.sh
    cd ..
done
