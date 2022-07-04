#!/bin/bash

rawdataset="/pubhome/qcxia02/dataset/TorsionNet500/TorsionNet500_mm_opt_geometries"

# for i in `seq 1 500`; do
for i in 5; do
    prefix="fragment_${i}_"
    count=$i

    python makedataset.py --prefix $prefix --count $count --rawdataset $rawdataset >> log
done