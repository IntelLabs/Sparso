#!/bin/bash

OUT_DIR=$1

for i in `cat pcg.lst`; do
bsub -C 0 -R "select[hsw] span[ptile=1]" -W 360 -q workq -J pcg ./pcg_cluster_child.sh $i $OUT_DIR
done
