#!/bin/bash

git submodule update --init

cp ./scripts/post-checkout .git/hooks 
cp ./scripts/post-rewrite .git/hooks 

./scripts/update-SpMP.sh
./scripts/build-julia.sh

