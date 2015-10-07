#!/bin/bash

for i in `cat pagerank.lst`; do
  echo !!!!$i
  echo %%%%auto-without-reordering
  julia pagerank-without-reordering.jl ~/matrices/$i.mtx
  echo
  echo %%%%auto
  julia pagerank.jl ~/matrices/$i.mtx
  echo
  echo
done
