#!/usr/bin/env bash

# Simple two-point plots
. PlotCommon.sh
PCMakeScriptOutDir
echo "$Ensemble: Three point - time reversed"
shopt -s nullglob

for PGroup in p2 #''
do
  for Spec in l s
  do (
    Dir=${Spec}${PGroup}
    InDir=$PlotData/corr/3pt_$Dir
    mkdir -p $Dir; cd $Dir
    for f in $InDir/quark_{l,s}_h*_g*_dt_*_p2_0_ps2_*_g*_g*.fold.*.txt; do
      if Split3ptFile $f; then
        [ $Gamma = "gXYZ" ] && Negate="y" || unset Negate
        Save=trev_${qSnk}_${qSrc}_${Gamma}_dt_${DeltaT}_p2_${pMax}_${opSnk}_${opSrc}
        if ! [ -f $Save.pdf ]
        then
          echo $Save
          if ! cmptrev.gp $f \
            $InDir/quark_${qSrc}_${qSnk}_${Gamma}_dt_${DeltaT}_p2_${ps2}_ps2_0_${opSrc}_${opSnk}.${f#*.} \
            $Save \
            "Time reverse check ${MSrcHuman}⟶${MSnkHuman}, n^2=${pMax}, ${GammaHuman} (${opHuman})" \
            "${MSrcHuman}⟶${MSnkHuman}" "${MSrcHuman}⟵${MSnkHuman}" "$Negate"
          then
            echo "Error $?"
          fi
        fi
      fi
    done
  ) done
done
