#!/usr/bin/env bash

# Simple two-point plots
. PlotCommon.sh
PCMakeScriptOutDir
echo "$Ensemble: Correlators two point"

# Allows selected plots to be turned off to make the script faster
DoCorr=$((1 - 0${corr+1})) # Plot the raw correlators, real and imaginary
DoRev=$((1 - 0${rev+1})) # compare real part of correlators with time-reversed self

Momenta=(0_0_0 1_0_0 1_1_0 1_1_1 2_0_0 2_1_0 2_1_1)

(( DoCorr )) && mkdir -p Corr

for (( pSq=0; pSq <= MaxPSq; ++pSq ))
do
  echo p ${Momenta[pSq]}
  # Temporarily, M1 zero momentum comes from special directory
  if [ "$Ensemble" == M1 ] && (( pSq == 0 ))
  then
    InDir=$PlotData/bootstrapf
  else
    InDir=$PlotData/bootstrap
  fi
  InDir=$InDir/2pt

  for f in $InDir/*_p_${Momenta[pSq]}_g*_g*.bootstrap.*.txt
  #for f in $InDir/s_l_p_1_?_?_g5P_g5P.bootstrap.*.txt
  do
    Filename=${f##*/}
    Extra=(${Filename//./ })
    base=${Extra[0]}
    ext=${Extra[-1]}
    random=${Extra[-2]}
    parts=(${base//_/ })
    if [ ${#Extra[@]} = 4 ] && [ ${#parts[@]} = 8 ] && [ "${parts[2]}" = p ]; then
      meson=${parts[0]}_${parts[1]}
      PCQNameNoMass "${parts[1]}" qOne
      PCQNameNoMass "${parts[0]}" qTwo
      MesonName="${qTwo}-${qOne}"
      snk=${parts[-2]}
      src=${parts[-1]}
      PCPointWall $snk $src opHuman
      # Real and imaginary component of each correlator
      if (( DoCorr )); then
        mkdir -p Corr/$meson
        title="$MesonName p(${parts[3]}, ${parts[4]}, ${parts[5]}) $opHuman (sink-source)"
        Legend="Real Imaginary"
        save=Corr/$meson/$base
        if ! [ -f $save.pdf ]
        then
          save=$save key='top right' \
            fields='corr corr_im' legend="$Legend" \
            savelabel= title="$title" \
            xlabel='t/a' ylabel="C(t/a)" \
            plot.sh $f
        fi
      fi
      # compare real part of correlators with time-reversed self
      if (( DoRev )) && (( pSq <= 1 )); then
        mkdir -p Rev/$meson
        unset Neg; [ "${snk:0:-1}" != "${src:0:-1}" ] && Neg='-'
        unset Final; [ "$snk" == g5W ] && [ "$Src" == gT5P ] && Final=29
        save=Rev/$meson/$base
        if ! [ -f $save.pdf ]
        then
          if ! cmptrev.gp $f '' $save '' '' '' "$Neg" '' $Final
          then
            echo "Error $? making $save.pdf from $f"
          fi
        fi
      fi
    fi
  done
done
