#!/usr/bin/env bash

export Ensemble=${Ensemble:-C2}
. FitZV.sh

#set -x

# Optional y range for plot of each ZV per Action,Spectator
aYRange[h${Heavy},h${Heavy}]=1.075:1.125
aYRange[h${Heavy},s]=0.97:1.1
aYRange[h${Heavy},l]=0.9:1.3 # Full range
#aYRange[h${Heavy},l]=0.97:1.05 # Zoom
aYRange[l,h${Heavy}]=0.65:0.9
aYRange[l,s]=0.71:0.75
aYRange[l,l]=0.71:0.75
aYRange[s,h${Heavy}]=0.68:0.76
aYRange[s,l]=0.71:0.78

if [ -v DoEFit ]; then
  [ -e $ZVEFit ] && rm $ZVEFit
  ti=8 tf=24 ti2=13 e2=1 PlotTI='8 8' PlotTF='28 28' FitZVEnergy h${Heavy} s 1.095:1.135
  ti=5 tf=20 ti2=7  e2=1 PlotTI='4 4' PlotTF='28 28' FitZVEnergy s l 0.315:0.36
  ti=8 tf=21 ti2=7 PlotTI='7 7' PlotTF='28 28' FitZVEnergy h${Heavy} l 1.05:1.12
  ti=12 tf=28 PlotTI='7 7' PlotTF='30 30' FitZVEnergy h${Heavy} h${Heavy} 1.64:1.665
  ti=6 tf=20 ti2=9 tf2=13 e2=1 PlotTI='4 4' PlotTF='28 28' FitZVEnergy s s 0.38:0.415
  ti=11 tf=20 ti2=6 e=1 PlotTI='4 4' PlotTF='28 28' FitZVEnergy l l 0.225:0.265
fi

if [ -v DoZV ]; then MakeZV; fi

if [ -v DoPlot ]; then PlotZV; fi

# Original ZV fits from constant
function FitConst()
{
(
  export UnCorr=
  ti=8 tf=16 FitZV h${Heavy} 24 s # Preferred
  ti=5 tf=12 FitZV l 16 s # Preferred
  unset ZVFit # End of preferred fit choices
  ti=5 tf=7 FitZV h${Heavy} 12 l
  ti=7 tf=9 FitZV h${Heavy} 16 l
  ti=8 tf=12 FitZV h${Heavy} 20 l # AltZV
  ti=8 tf=16 FitZV h${Heavy} 24 l
  ti=8 tf=12 FitZV h${Heavy} 20 s
  ti=9 tf=19 FitZV h${Heavy} 28 s
  ti=8 tf=24 FitZV h${Heavy} 32 s

  ti=7 tf=18 FitZV l 24 s
  ti=5 tf=13 FitZV l 16 l
  ti=7 tf=19 FitZV l 24 l
)
}

function FitJan24Light()
{
(
  Q=l
  Spec=s
  yrange=0.72:0.75
  e=3 ti='6 8 10 12' tf='10 12 14 16' FitZVNew 16 20 24 28 # Preferred
  unset ZVFit # End of preferred fit choices
  e=3 ti='6 8 10 12' tf='10 12 14 16' UnCorr= FitZVNew 16 20 24 28
  (
  e=3
  for((i=0;i<2;++i)); do
    ti='6 8' tf='10 12' FitZVNew 16 20
    #ti='5 5 5' tf='11 15 19' FitZVNew 16 20 24 # Offset
    ti='8 10' tf='12 14' FitZVNew 20 24
    #ti='5 5' tf='15 19' FitZVNew 20 24 # offset
    # dT=12 is just too early
    #ti='4 5 5' tf='8 11 15' FitZVNew 12 16 20 # Shite
    #ti='5 6 8' tf='7 11 13' FitZVNew 12 16 20 # Also shite
    UnCorr=
  done
  )
  e=2
  #ti=3 tf=9 Stat=0.02 FitZVNew 12
  ti=4 tf=8 FitZVNew 12
  ti=5 tf=11 FitZVNew 16
  ti=5 tf=15 FitZVNew 20
  ti=5 tf=19 FitZVNew 24
)
}

function FitJan24Heavy()
{
(
  Q=h${Heavy}
  Spec=s
  yrange=0.98:1.06
  e=3 ti='8 9' tf='12 15' FitZVNew 20 24 # Preferred
  unset ZVFit # End of preferred fit choices
  ti=8 tf=16 FitZVNew 24 # Best of single correlators
  ti=8 tf=12 FitZVNew 20 # 2nd-best 1-corr (smaller errors)
  ti=9 tf=11 FitZVNew 20
  ti=9 tf=15 FitZVNew 24
  ti=10 tf=14 FitZVNew 24
  export UnCorr=
  ti='6 8' tf='10 12' FitZVNew 16 20
  ti='6 8 9' tf='10 12 14' FitZVNew 16 20 24
  ti='8 9' tf='12 14' FitZVNew 20 24
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
    ZVLink=ZVAltZV.txt # Alternate ZV choices are actually the original choices on C2
  done
  )
fi
