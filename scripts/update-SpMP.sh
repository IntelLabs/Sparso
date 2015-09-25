#!/bin/bash

HTTP_PROXY="" git submodule update deps/SpMP

cd deps/SpMP

git checkout master

# update
git pull

cd ../..

