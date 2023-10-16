#!/usr/bin/env bash

# Walk all the ratios looking for error at maximum
ratio=${ratio:-ratio}

# Plot ratios
. PlotCommon.sh
#PCMakeScriptOutDir
OptionASCII=

#set -x

echo "$Ensemble: Find correlator maximum (and error)"
shopt -s nullglob

for d in $PlotData/$ratio/*
do
    Dir=${d##*/}
    if [ "$Dir" != log ]
    then
      # Get spectator and how momenta organised
      unset PGroup
      Spec=$Dir
      if [ "${Spec: -2}" == "p2" ]
      then
        PGroup=${Spec: -2}
        Spec=${Spec:0:-2}
      fi
      # Walk file list
      for f in $d/R*.fold.*.txt
      do
        if Split3ptFile $f
        then
          SaveName=${f##*/}
          SaveName=${SaveName/.fold./.errmax.}
          echo ${Ratio:1} "\"$MSnkHuman\"" "\"$MSrcHuman\"" $Gamma $DeltaT $pMax $opSnk $opSrc \
          $(gawk 'BEGIN {Max=-9E99}; /^[-0-9]/ { if ($4 > Max) { t=$1; Low=$3; Max=$4; High=$5 }}; END { print t, (High - Low)/(2*Max), Low, Max, High }' $f)
        fi
      done
    fi
done
