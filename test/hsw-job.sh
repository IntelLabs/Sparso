#!/bin/bash

N=1
if [ "$#" -ge "1" ] ; then
  N=$1
fi

bsub -n $N -C 0 -R "select[hsw] span[ptile=1]" -W 720 -q workq -o /panfs/users/hrong1/SparseAccelerator/test/p.out -e /panfs/users/hrong1/SparseAccelerator/test/p.err ./perf.sh

