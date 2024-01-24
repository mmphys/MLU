#!/usr/bin/env bash

export Ensemble=${Ensemble:-M2}
. FitZV.sh

#set -x

# Optional y range for plot of each ZV per Action,Spectator
#aYRange[h${Heavy},h${Heavy}]=1.06:1.13
#aYRange[h${Heavy},s]=0.81:0.845
#aYRange[h${Heavy},l]=0.9:1.5
#aYRange[l,h${Heavy}]=0.6:1
#aYRange[l,s]=0.69:0.73
#aYRange[l,l]=0.69:0.73
#aYRange[s,h${Heavy}]=0.68:0.76
#aYRange[s,l]=0.695:0.735

if [ -v DoEFit ]; then
  [ -e $ZVEFit ] && rm $ZVEFit
  PlotTF='29 29'
  ti=13 tf=28 ti2=18 PlotTI='10 10' FitZVEnergy h${Heavy} h${Heavy} 1.235:1.253
  ti=10 tf=26 ti2=19 e2=1 PlotTI='7 11' FitZVEnergy h${Heavy} s 0.815:0.845
  ti=9 tf=21 ti2=18 tf2=23 e2=1 PlotTI='7 7' FitZVEnergy h${Heavy} l 0.77:0.84
  #ti=9 tf=21 ti2=12 tf2=21 e2=1 PlotTI='7 7' FitZVEnergy h${Heavy} l 0.77:0.84
  ti=7  tf=20 ti2=12 e2=1 PlotTI='6 6' FitZVEnergy s s .28:.34
  ti=6 tf=20 ti2=8 e2=1 PlotTI='5 5' FitZVEnergy s l 0.22:0.28
  ti=6 tf=20      e2=1 PlotTI='5 5' FitZVEnergy l l .13:.19
  #ti=9 tf=20 ti2=6 e=1 PlotTI='5 5' FitZVEnergy l l .13:.19
fi

if [ -v DoZV ]; then MakeZV; fi

if [ -v DoPlot ]; then PlotZV; fi

# Original ZV fits from constant
function FitConst()
{
(
  [ -e $ZVFit ] && unset ZVFit #Don't overwrite fit selection if it already exists
  export UnCorr=
  ti=7 tf=9 FitZV h${Heavy} 16 s
  ti=8 tf=12 FitZV h${Heavy} 20 s
  ti=8 tf=14 FitZV h${Heavy} 24 s # Preferred
  ti=8 tf=20 FitZV h${Heavy} 28 s
  ti=7 tf=9 FitZV h${Heavy} 16 l
  ti=7 tf=12 FitZV h${Heavy} 20 l # AltZV
  ti=5 tf=18 FitZV h${Heavy} 24 l

  ti=4 tf=12 FitZV l 16 s
  ti=4 tf=16 FitZV l 20 s # Preferred
  ti=4 tf=20 FitZV l 24 s
  ti=8 tf=20 FitZV l 28 s
  ti=4 tf=12 FitZV l 16 l
  ti=4 tf=16 FitZV l 20 l
  ti=4 tf=20 FitZV l 24 l
  ti=4 tf=24 FitZV l 28 l
)
}

function FitJan24Light()
{
(
  Q=l
  Spec=s
  yrange=0.74:0.775
  e=3 ti='6 8' tf='10 12' FitZVNew 16 20 # Preferred
  unset ZVFit # End of preferred fit choices
  e=3 ti='6 8 10 12' tf='10 12 14 16' FitZVNew 16 20 24 32 # Preferred
)
}

function FitJan24Heavy()
{
(
  Q=h${Heavy}
  Spec=s
  yrange=0.98:1.04
  e=3 ti='8 10 12 14' tf='12 14 16 18' FitZVNew 20 24 28 32 # Preferred
  unset ZVFit # End of preferred fit choices
  e=3 ti='8 10 12' tf='12 14 16' FitZVNew 20 24 28
  e=3 ti='8 10' tf='12 14' FitZVNew 20 24
  #e=3 ti='6 8 10' tf='10 12 14' FitZVNew 16 20 24
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
  Q=h${Heavy}
  Spec=s
  yrange=0.98:1.04
  # Test thinning
  dT=28
  export ti=8
  export tf=$((dT-ti))
  for Thin in {1..12}; do
    Thin=$Thin ExtraName=thin_$Thin FitZV h${Heavy} $dT s
  done
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
