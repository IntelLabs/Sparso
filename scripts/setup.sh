#!/bin/bash

git submodule update --init

echo "Setting up hooks"
cp ./scripts/post-checkout .git/hooks 
cp ./scripts/post-rewrite .git/hooks 

echo "Building SpMP lib"
./scripts/update-SpMP.sh
#/opt/intel/bin/compilervars.sh intel64
cd ./lib
make clean
make -j DBG=yes
cd ..

echo "Patching CompilerTools"
cd deps/CompilerTools/
patch -p 1 < ../../scripts/CompilerTools.patch
cd - 

