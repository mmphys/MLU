#!/usr/bin/env bash

for opsnk in g5P; do
for opsrc in g5P g5W; do

# Individual files by Q^2
outDir=byq2
mkdir -p $outDir
for hsrc in {1..3}; do
for hsnk in `seq 0 $((hsrc - 1))`; do
save=$outDir/h${hsnk}_${opsnk}_h${hsrc}_${opsrc} title="h${hsnk},${opsnk} <- h${hsrc},${opsrc}" xargs='12 14 16 20' key='bottom center' x='column(1)/word(xargs,File)' offset=0.001 ti=-0.01 tf=1.01 fields='corr' plot.sh 3pt_s/R2_h${hsnk}_p_0_0_0_h${hsrc}_p_0_0_0_dt_{12,14,16,20}_${opsnk}_${opsrc}.fold.1835672416.txt
done
done

# Individual plots by DeltaT
outDir=bydt
mkdir -p $outDir
for dt in 12 14 16 20; do
save=$outDir/h0_dt_${dt}_${opsnk}_${opsrc} title="DeltaT=${dt}, h0 ${opsnk} <- h_n ${opsrc}" ti=0 tf=${dt} yrange='*:*' fields='corr' plot.sh 3pt_s/R2_h0_p_0_0_0_h*_p_0_0_0_dt_${dt}_${opsnk}_${opsrc}.fold.1835672416.txt
done

done
done
