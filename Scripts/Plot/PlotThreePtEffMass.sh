#!/usr/bin/env bash

# Simple two-point plots
. PlotCommon.sh
PCMakeScriptOutDir
echo "$Ensemble: Three point - effective mass"
shopt -s nullglob

PlotFields=(exp cosh)
PlotFieldTitle=('log' 'cosh')

for field in {0..1}
do
  if [ "${PlotFields[field]}" == exp ]; then
    negate="-"
    NegOffset=1
  else
    unset negate
    NegOffset=0
  fi
  for PGroup in p2 #''
  do
    for Spec in l s
    do (
      Dir=${Spec}${PGroup}
      InDir=$PlotData/corr/3pt_$Dir
      Dir=${PlotFieldTitle[field]}/$Dir
      echo Dir="\"$Dir\""
      mkdir -p $Dir; cd $Dir
      for f in $InDir/quark_{l,s}_h*_g*_dt_*_p2_0_ps2_*_g*_g*.fold.*.txt; do
        if Split3ptFile $f; then
          save=${Ratio}_${qSnk}_${qSrc}_${Gamma}_dt_${DeltaT}_p2_${p2}_ps2_${ps2}_${opSnk}_${opSrc}_${PlotFieldTitle[field]} \
            legend="'${MSrcHuman}⟶${MSnkHuman}' '${negate}${MSrcHuman}⟵${MSnkHuman}'" \
            fields=${PlotFields[field]} ti=1.5 tf=$((DeltaT-2)).5 savelabel= \
            title="${PlotFieldTitle[field]} effective mass ${MSrcHuman}⟶${MSnkHuman}, n^2=${pMax}, ${GammaHuman} (${opHuman})" \
            xlabel=t/a ylabel="${PlotFieldTitle[field]} effective mass" \
            negate="+ $negate" \
            xrevt="0 $((DeltaT+NegOffset))" \
            plot.sh $f \
            $InDir/quark_${qSrc}_${qSnk}_${Gamma}_dt_${DeltaT}_p2_${ps2}_ps2_0_${opSrc}_${opSnk}.${f#*.}
        fi
      done
    ) done
  done
done
