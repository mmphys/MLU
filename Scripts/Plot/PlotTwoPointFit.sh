#!/usr/bin/env bash

# Simple two-point plots
. PlotCommon.sh
PCMakeScriptOutDir
echo "$Ensemble: Two point fits"

InDir=$PlotData/fit/2ptp2/Summary

# Optional environment variables
# nonew:  Set to anything to disable new plots
# noold:  Set to anything to disable old plots
DoNew=$((1-0${nonew+1}))
DoOld=$((1-0${noold+1}))

declare -A WhereClauses
WhereClauses[eq]='column("ti") == column("ti1")'
WhereClauses[lt]='column("ti") <  column("ti1")'
WhereClauses[gt]='column("ti") >  column("ti1")'
extralbl='column("ti") == column("ti1") ? "" : "(".stringcolumn("ti1").")"'

for f in $InDir/*.params_sort.*.txt
do
  Dir=${f%/*}
  Filename=${f##*/}
  Ext=${Filename#*.}
  Prefix=${Filename%%.*}
  PrefixParts=(${Prefix//_/ })
  if [ ${#PrefixParts[@]} == 4 ] && [ "${PrefixParts[2]}" == "p2" ]
  then
    echo $Prefix
    NumRows=$(awk '/^[-0-9]/ {z++}; END {print z}' < $f)
    if (( NumRows && DoOld ))
    then
      Cmd="s1plotparams.sh ${f/params_sort/params}"
      (( NumRows > 75 )) && Cmd="label= $Cmd"
      eval $Cmd
    fi
    if (( NumRows && DoNew ))
    then
      for series in tf ti
      do
        for WC in "${!WhereClauses[@]}"
        do
          #sort='-V -k 5r -k 1'
          series=$series save="${series}_$WC" where="${WhereClauses[$WC]}" \
            extralbl="$extralbl" sort=$sort PlotFitRes.sh $f
        done
      done
    fi
  fi
done
