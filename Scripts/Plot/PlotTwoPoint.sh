#!/usr/bin/env bash

# Simple two-point plots
. PlotCommon.sh
PCMakeScriptOutDir
echo "$Ensemble: Two point"

# Allows selected plots to be turned off to make the script faster
bypw="${bypw:-1}" # plots grouped by source / sink (point - wall)
bymom="${bymom:-1}" # plots grouped by momentum
relerr="${relerr:-1}" # Do relative error
masslog="${masslog:-1}" # Log effective mass
masscosh="${masscosh:-1}" # Cosh effective mass
if [ "$mass" = "0" ]; then
  masslog=0
  masscosh=0
fi

domass=("$masslog" "$masscosh")

InDir=$PlotData/corr/2ptp2

PlotFields=(exp cosh)
PlotFieldTitle=('log ' 'cosh ')

OutSubDirs=(All Zoom)
OutTI=(4.5 7.5)
#OutTI=(4.5 4.5)
OutTIRelErr=(0.5 4.5)
OutTF=($((THalf-1)).5 $((THalf*6/10-1)).5)
#OutTF=($((THalf-1)).5 23.5)
WhichSinks=('g5P g5W' g5P)

for Dir in {0..1}
do (
  mkdir -p ${OutSubDirs[Dir]}; cd ${OutSubDirs[Dir]}
  for f in $InDir/*_p2_0_g5P_g5P.fold.1835672416.txt
  do
    base=${f##*/}
    ext=${base#*.}
    base=${base%%.*}
    parts=(${base//_/ })
    #echo $base
    if [ ${#parts[@]} = 6 ]; then
      OptionSimpleTitle= # Comment this out for more info in title
      meson=${parts[0]}_${parts[1]}
      PCQNameNoMass "${parts[1]}" qOne
      PCQNameNoMass "${parts[0]}" qTwo
      OptionNoMass=
      GetMeson MesonName "${qTwo}" "${qOne}"
      # Plots by source and sink with all momenta on same chart
      if (( bypw )); then
        for snk in g5P g5W; do
          for src in g5P g5W; do
            PCPointWall $snk $src opHuman
            title="$MesonName"
            if ! [ -v OptionSimpleTitle ]; then title="$title $opHuman (sink-source)"; fi
            unset FileNames
            unset Legend
            pSq=$(( MaxPSq + 1 ))
            while (( pSq-- ))
            do
              FileNames="$FileNames $InDir/${meson}_p2_${pSq}_${snk}_${src}.$ext"
              Legend="$Legend n^2=$pSq"
            done
            [ $snk == $src ] && key='top right' || key='bottom right'
            for field in {0..1}; do
              if (( domass[field] )); then
                [ -v OptionSimpleTitle ] && ylabel="a E_{eff}" || ylabel="${PlotFieldTitle[field]}effective mass"
                save=${meson}_${snk}_${src}_${PlotFields[field]} \
                  fields=${PlotFields[field]} ti=${OutTI[Dir]} tf=${OutTF[Dir]} key="$key" legend="$Legend" \
                  savelabel= title="$title" \
                  xlabel='t/a' ylabel="$ylabel" legenh= \
                  plot.sh $FileNames
              fi
            done
            [ $snk == g5W ] && key='top left' || key='bottom right'
            if (( relerr )); then
              save=${meson}_${snk}_${src}_relerr \
                fields=corr ti=${OutTIRelErr[Dir]} tf=${OutTF[Dir]} key="$key" legend="$Legend" \
                rel=2 log=1 savelabel= title="$title" \
                xlabel='t/a' ylabel="Relative error" \
                plot.sh $FileNames
            fi
          done
        done
      fi
      # Plots by momenta with all sources and sinks on same chart
      if (( bymom )); then
        pSq=$(( MaxPSq + 1 ))
        while (( pSq-- ))
        do
          title="$MesonName n^2=$pSq (sink-source)"
          (( PSq <= MaxPSq / 2 )) && key='top right' || key='bottom right'
          unset FileNames
          unset Legend
          for snk in ${WhichSinks[Dir]}; do
            for src in g5P g5W; do
              PCPointWall $snk $src opHuman
              FileNames="$FileNames $InDir/${meson}_p2_${pSq}_${snk}_${src}.$ext"
              Legend="$Legend $opHuman"
            done
          done
          for field in {0..1}; do
            if (( domass[field] )); then
              save=${meson}_p2_${pSq}_${PlotFields[field]} \
                fields=${PlotFields[field]} ti=${OutTI[Dir]} tf=${OutTF[Dir]} key="$key" legend="$Legend" \
                savelabel= title="$title" \
                xlabel='t/a' ylabel="${PlotFieldTitle[field]}effective mass" \
                plot.sh $FileNames
            fi
          done
          if (( relerr )); then
            save=${meson}_p2_${pSq}_relerr \
              fields=corr ti=${OutTIRelErr[Dir]} tf=${OutTF[Dir]} key='bottom right' legend="$Legend" \
              rel=2 log=1 savelabel= title="$title" \
              xlabel='t/a' ylabel="Relative error" \
              plot.sh $FileNames
          fi
        done
      fi
    fi
  done
) done
