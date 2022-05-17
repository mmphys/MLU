#!/usr/bin/env bash

# Optional environment variables
# dir: subdirectory, e.g."frozen"

# Plot ratios
. PlotCommon.sh
PCMakeScriptOutDir

#set -x

# Computed from input
OutSub=${dir:-final}
mkdir -p $OutSub; cd $OutSub

echo "$Ensemble: Three point $OutSub ratios"
shopt -s nullglob

for RatioNum in {1..3}
do
  [ RatioNum == 3 ] && Log=1 || Log=0
  for PGroup in p2 #''
  do
    for Spec in l s
    do (
      Dir=${Spec}${PGroup}
      InDir=$PlotData/ratio$dir/3pt_$Dir
      Dir=R$RatioNum/$Dir
      echo Dir="\"$Dir\""
      mkdir -p $Dir; cd $Dir
      for f in $InDir/R${RatioNum}_*_g*P_g*P.fold.*.txt; do
        if Split3ptFile $f $Spec; then
          FilePrefix=${Ratio}_${qSnk}_${qSrc}_${Gamma}_dt_${DeltaT}_p2_${p2}
          unset FileNames
          unset Legend
          for snk in ${opSnk} ${opSnk%P}W; do
            for src in ${opSnk} ${opSnk%P}W; do
              PCPointWall $snk $src opHuman
              FileNames="$FileNames $InDir/${FilePrefix}_${snk}_${src}.$Ext"
              Legend="$Legend $opHuman"
            done
          done
          save=${FilePrefix} \
            legend="$Legend" \
            fields=corr ti=0.5 tf=$((DeltaT-1)).5 savelabel= \
            title="${Ratio} ${MSrcHuman}⟶${MSnkHuman}, n^2=${p2}, ${GammaHuman} (ΔT=${DeltaT})" \
            xlabel=t/a ylabel="${Ratio}" \
            plot.sh $FileNames
        fi
      done
    )
    done
  done
done
