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

echo "* Setting up BinDeps"
cd deps
julia build.jl
cd ..

echo "* Setting up submoduls"
git submodule update --init
cd deps/CompilerTools
git checkout 99656518b2396b8e6f92cb2cd6d23cffb5e3f8ad
cd ../..

echo "* Setting up git hooks"
cp ./scripts/post-checkout .git/hooks 
cp ./scripts/post-rewrite .git/hooks 

if [ -f ~/.julia/lib/v0.4/CompilerTools.ji ]; then
    echo "* Removing precompiled CompilerTools"
    rm ~/.julia/lib/v0.4/CompilerTools.ji
fi

echo "* Patching CompilerTools"
cd deps/CompilerTools/
patch -p 1 < ../../scripts/CompilerTools.patch
cd - 

echo "* Building SpMP lib"
. ./scripts/update-SpMP.sh
cd ./lib
make clean
make -j 
cd ..

echo "* Downloading extra matrices"
cd test
cd matrices
wget ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc2/bcsstk14.mtx.gz
gunzip bcsstk14.mtx.gz
wget https://raw.githubusercontent.com/exmatex/CoSP2/master/data/hmatrix.512.mtx
wget https://raw.githubusercontent.com/exmatex/CoSP2/master/data/hmatrix.1024.mtx

cd lp
for i in `cat ../../correctness/ipm.lst`; do
  wget http://www.cise.ufl.edu/research/sparse/MM/LPnetlib/${i}.tar.gz
  tar xzvf ${i}.tar.gz
  rm ${i}.tar.gz
done
cd ..

# TODO: download matrices of LBFGS, PageRank, and PCG

cd ..
cd ..
