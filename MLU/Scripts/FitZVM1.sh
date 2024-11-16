#!/usr/bin/env bash

export Ensemble=${Ensemble:-M1}
. FitZV.sh

#set -x

# Optional y range for plot of each ZV per Action,Spectator
#aYRange[h${Heavy},h${Heavy}]=1.06:1.13
#aYRange[h${Heavy},s]=0.96:1.02
#aYRange[h${Heavy},l]=0.9:1.5
#aYRange[l,h${Heavy}]=0.6:1
#aYRange[l,s]=0.72:0.84
#aYRange[l,l]=0.69:0.73
#aYRange[s,h${Heavy}]=0.68:0.76
#aYRange[s,l]=0.695:0.735

if [ -v DoEFit ]; then
  [ -e $ZVEFit ] && rm $ZVEFit
  ti=14 tf=28 ti2=19 FitZVEnergy h${Heavy} h${Heavy} 1.24:1.25
  #ti=9 tf=28 ti2=15 tf2=24 e2=1 PlotTI='7 11' FitZVEnergy h${Heavy} s 0.815:0.84
  ti=10 tf=28 ti2=15 tf2=28 e2=1 PlotTI='7 11' FitZVEnergy h${Heavy} s 0.815:0.84
  ti=9 tf=28 ti2=14 tf2=25 e2=1 PlotTI='7 11' FitZVEnergy h${Heavy} l 0.78:0.85
  ti=6  tf=22 ti2=12 e2=1 FitZVEnergy s s .28:.34
  ti=6 tf=18 ti2=11 e2=1 PlotTI='5 5' FitZVEnergy s l 0.21:0.26
  ti=9  tf=20 ti2=6  e=1  FitZVEnergy l l .1:.14
fi

if [ -v DoZV ]; then MakeZV; fi

if [ -v DoPlot ]; then PlotZV; fi

# Original ZV fits from constant
function FitConst()
{
(
  export UnCorr=
  ti=8 tf=15 FitZV h${Heavy} 24 s # Preferred
  ti=6 tf=11 FitZV l 16 s # Preferred ZV_l
  unset ZVFit # End of preferred fit choices
  ti=7 tf=9 FitZV h${Heavy} 16 s
  ti=8 tf=12 FitZV h${Heavy} 20 s
  ti=8 tf=20 FitZV h${Heavy} 28 s
  ti=6 tf=9 FitZV h${Heavy} 16 l
  ti=8 tf=12 FitZV h${Heavy} 20 l # AltZV
  ti=8 tf=15 FitZV h${Heavy} 24 l
  ti=8 tf=20 FitZV h${Heavy} 28 l

  ti=5 tf=15 FitZV l 20 s
  ti=5 tf=19 FitZV l 24 s
  ti=5 tf=23 FitZV l 28 s
  ti=6 tf=11 FitZV l 16 l
  ti=5 tf=15 FitZV l 20 l
  ti=5 tf=19 FitZV l 24 l
  ti=5 tf=23 FitZV l 28 l
)
}

function FitJan24Light()
{
(
  Q=l
  Spec=s
  yrange=0.73:0.76
  e=3 ti='6 8' tf='10 12' FitZVNew 16 20 # Preferred
  unset ZVFit # End of preferred fit choices
)
}

function FitJan24Heavy()
{
(
  Q=h${Heavy}
  Spec=s
  yrange=0.97:1
  e=3 ti='8 10' tf='12 14' FitZVNew 20 24 # Preferred
  unset ZVFit # End of preferred fit choices
  export UnCorr=
)
}

# Jan 2024 ZV Fits from model including excited-states
function FitJan24()
{
  FitJan24Light
  FitJan24Heavy
}

function FitTest()
{
(
  Q=l
  Spec=s
  yrange=0.72:0.75
  e=3 ti='6 8 10 12' tf='10 12 14 16' FitZVNew 16 20 24 28 # Preferred
)
}

############################################################

# Choose which ZV to perform. Defaults to DefaultSeries (see FitZV.sh)

############################################################

if [ -v series ] || [ -v DoFit ]; then
  for Series in ${series:-Const Jan24}
  do
    Series=${Series@L}
    Series=${Series@u}
    case $Series in
      Const | Jan24 )
        ZVFit=ZV$Series.txt
        [ -e $ZVFit ] && rm $ZVFit #Overwrite fit selection if it already exists
        eval Fit$Series;;
      Test ) unset ZVFit; FitTest;;
      *) echo "Ignoring unknown series '$Series'";;
    esac
  done
fi

############################################################

# Activate the default ZV series

############################################################

ZVFit=ZV$DefaultSeries.txt
if ! [ -e "$ZVFit" ]; then
  echo "Can't link $ZVLink -> '$ZVFit' (doesn't exist)"
else
  (
  for ((i=0;i<2;++i)); do
    [ -e "$ZVLink" ] || [ -h "$ZVLink" ] && rm "$ZVLink"
    ln -s "$ZVFit" "$ZVLink"
    ZVLink=ZVAltZV.txt # Alternate ZV choices are actually the original choices on F1M
  done
  )
fi
