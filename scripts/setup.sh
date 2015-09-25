#!/bin/bash

git submodule update --init

echo "Setting up hooks"
cp ./scripts/post-checkout .git/hooks 
cp ./scripts/post-rewrite .git/hooks 

#echo "Setting up julia"
#wget https://julialang.s3.amazonaws.com/bin/linux/x64/0.4/julia-0.4.0-rc1-linux-x86_64.tar.gz

echo "Building SpMP lib"
./scripts/update-SpMP.sh
/opt/intel/bin/compilervars.sh intel64
cd ./lib
make clean
make -j DBG=yes
cd ..

#echo "Building julia"
#./scripts/build-julia.sh
