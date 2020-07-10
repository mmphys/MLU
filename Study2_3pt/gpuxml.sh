#!/bin/bash

for MyType in GFPW Z2; do
for (( t = 0 ; t < 64 ; t += 8 )); do
  jobxml=GPU_${MyType}_${t}.xml
  sed -e "s|@type@|$MyType|g" -e "s|@tstart@|$t|g" -e "s|@tend@|$((t+8))|g" xml3pt.xml > $jobxml
  xml3pt $jobxml
  rm $jobxml
done
done
