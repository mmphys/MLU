#!/usr/bin/env bash

export Ensemble=${Ensemble:-F1M}
. FitZV.sh

#set -x

# Optional y range for plot of each ZV per Action,Spectator
aYRange[h${Heavy},h${Heavy}]=0.995:1.02
aYRange[h${Heavy},s]=0.95:1.01
aYRange[h${Heavy},l]=0.96:1.1
aYRange[l,h${Heavy}]=0.75:0.86
aYRange[l,s]=0.75:0.78
aYRange[l,l]=0.745:0.775
aYRange[s,h${Heavy}]=0.74:0.782
aYRange[s,l]=0.75:0.78

if [ -v DoEFit ]; then
  [ -e $ZVEFit ] && rm $ZVEFit
  #PlotTF='46 46'
  #ti=6 tf=29 ti2=5 e=3 PlotTI='5 5' FitZVEnergy h${Heavy} s 0.68:0.83 # Old (bad) choice
  ti=16 tf=40 tf2=29 FitZVEnergy h${Heavy} s 0.723:0.736
  #ti=10 tf=26 ti2=7 PlotTI='7 5' FitZVEnergy s l 0.18:0.22 # Old (bad) choice
  ti=8 tf=31 ti2=8 tf2=20 e2=1 PlotTI='7 5' FitZVEnergy s l 0.18:0.22
  #ti=8 tf=19 ti2=8 tf2=20 e2=1 PlotTI='7 5' FitZVEnergy s l 0.18:0.22
  ti=10 tf=27 PlotTI='7 7' FitZVEnergy h${Heavy} l 0.68:0.726
  ti=16 tf=37 ti2=29 tf2=44 e2=1 PlotTI='10 10' FitZVEnergy h${Heavy} h${Heavy} 1.098:1.105
  #ti=16 tf=37 PlotTI='10 10' FitZVEnergy h${Heavy} h${Heavy} 1.094:1.105
  ti=10 tf=25 ti2=6 PlotTI='9 6' FitZVEnergy s s 0.248:0.269
  ti=7 tf=24 ti2=7 tf2=23 e2=1 PlotTI='7 5' FitZVEnergy l l .07:0.1
  unset PlotTF
fi

if [ -v DoZV ]; then MakeZV; fi

if [ -v DoPlot ]; then PlotZV; fi

# Original ZV fits from constant
function FitConst()
{
(
  export UnCorr=
  ti=10 tf=14 FitZV h${Heavy} 24 s # Preferred
  ti=4 tf=16 FitZV l 20 s # Preferred
  unset ZVFit # End of preferred fit choices
  ti=10 tf=18 FitZV h${Heavy} 28 s
  ti=10 tf=22 FitZV h${Heavy} 32 s
  ti=7 tf=10 FitZV h${Heavy} 16 l
  ti=8 tf=13 FitZV h${Heavy} 20 l # AltZV
  ti=8 tf=18 FitZV h${Heavy} 24 l
  ti=8 tf=20 FitZV h${Heavy} 28 l
  ti=8 tf=26 FitZV h${Heavy} 32 l

  ti=4 tf=20 FitZV l 24 s
  ti=4 tf=24 FitZV l 28 s
  ti=4 tf=28 FitZV l 32 s
  ti=4 tf=16 FitZV l 20 l
  ti=4 tf=20 FitZV l 24 l
  ti=4 tf=24 FitZV l 28 l
  ti=4 tf=28 FitZV l 32 l
)
}

function FitJan24Light()
{
(
  Q=l
  Spec=s
  yrange=0.755:0.77
  e=3 ti='6 8 10 12' tf='10 12 14 16' FitZVNew 16 20 24 28 # Preferred
  unset ZVFit # End of preferred fit choices
)
}

function FitJan24Heavy()
{
(
  Q=h${Heavy}
  Spec=s
  yrange=0.97:1
  e=3 ti='8 9' tf='12 15' FitZVNew 20 24 # Preferred
  unset ZVFit # End of preferred fit choices
  e=3 ti='6 8 9' tf='10 12 15' FitZVNew 16 20 24 # Preferred
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
