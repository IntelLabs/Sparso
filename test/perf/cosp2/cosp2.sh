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

for i in hmatrix.512.mtx hmatrix.1024.mtx; do
  echo "#### INPUT: $i"

  input=../inputs/${i}
  if [[ ! -f $input ]]; then
    echo "SKIP: file '${input}' does't exist"
    echo
    continue
  fi
  echo

  #OMP_NUM_THREADS=1 julia cosp2.jl ../matrices/$i.mtx
  echo "---- BEGIN: Julia-as-is"
  julia cosp2.jl $input "julia"
  echo "---- END: Julia-as-is"
  echo

  echo "---- BEGIN: Baseline (Call-repl)"
  julia cosp2.jl $input "call-repl"
  echo "---- END: Baseline"
  echo

  echo "---- BEGIN: +Matrix-properties"
  julia cosp2.jl $input "context"
  echo "---- END: +Matrix-properties"
  echo
done
