#!/bin/bash
git submodule update

cd deps/julia
sed -n 's/1048576/16777216/' base/pcre.jl
make -j 4 OPENBLAS_TARGET_ARCH=NEHALEM debug
cd ../..
