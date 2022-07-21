#!/usr/bin/env bash

# Optional environment variables
# dir: subdirectory, e.g."frozen"

# Plot ratios
. PlotCommon.sh
PCMakeScriptOutDir

#set -x

# Computed from input
OutSub=${dir:-final}

# Optional environment variables
# nodt:  Set to anything to disable plots individual DeltaT all point-wall
# nopw:  Set to anything to disable plots individual point-wall all DeltaT
DoDeltaT=$((1-0${nodt+1}))
DoPW=$((1-0${nopw+1}))

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
      InDir=$PlotData/ratio$dir/$Dir #3pt_$Dir
      Dir=R$RatioNum/$Dir
      echo Dir="\"$Dir\""
      mkdir -p $Dir; cd $Dir
      # For each Delta T, show the various point/wall sources/sinks
      if (( DoDeltaT ))
      then
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
      fi
      # For each souece/sink (point/wall), show the various Delta T
      if (( DoPW ))
      then
        for f in $InDir/R${RatioNum}_*_dt_*_g*_g*.fold.*.txt; do
          if Split3ptFile $f $Spec; then
            FilePrefix=${Ratio}_${qSnk}_${qSrc}_${Gamma}
            FileSuffix=p2_${p2}_${opSnk}_${opSrc}
            PCPointWall $opSnk $opSrc opHuman
            unset FileNames
            unset Legend
            for DeltaT in ${EnsembleDeltaT[@]}; do
              FileNames="$FileNames $InDir/${FilePrefix}_dt_${DeltaT}_${FileSuffix}.$Ext"
              Legend="$Legend ΔT/a=${DeltaT}"
            done
            save=${FilePrefix}_${FileSuffix} \
              legend="$Legend" \
              fields=corr savelabel= \
              title="${Ratio} ${MSrcHuman}⟶${MSnkHuman}, n^2=${p2}, ${GammaHuman} (${opHuman})" \
              xlabel=(t-ΔT/2)/a ylabel="${Ratio}" offset=0 \
              x='((column(1)<2 || column(1)>word(FileDT,File)-2) ? NaN : column(1)-0.5*word(FileDT,File))' \
              plot.sh $FileNames
          fi
        done
      fi
    )
    done
  done
done
