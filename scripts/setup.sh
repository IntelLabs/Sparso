#!/bin/bash

echo "* Setting up submoduls"
git submodule update --init

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

echo "* Extracting matrices"
cd test
tar xf matrices.tar.gz
cd ..
