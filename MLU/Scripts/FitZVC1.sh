#!/usr/bin/env bash

export Ensemble=${Ensemble:-C1}
. FitZV.sh

#set -x

# Optional y range for plot of each ZV per Action,Spectator
aYRange[h${Heavy},h${Heavy}]=1.06:1.13
aYRange[h${Heavy},s]=0.95:1.06
aYRange[h${Heavy},l]=0.9:1.5
aYRange[l,h${Heavy}]=0.6:1
aYRange[l,s]=0.69:0.73
aYRange[l,l]=0.69:0.73
aYRange[s,h${Heavy}]=0.68:0.76
aYRange[s,l]=0.695:0.735

if [ -v DoEFit ]; then
  [ -e $ZVEFit ] && rm $ZVEFit
  ti=6 tf=27 ti2=14 e2=1 PlotTI='5 10' FitZVEnergy h${Heavy} s 1.095:1.12
  ti=6 tf=23 ti2=7  e2=1 PlotTI='5 5'  FitZVEnergy s l 0.30:0.325
  ti=8 tf=28 ti2=14 e2=1 PlotTI='7 7'  FitZVEnergy h${Heavy} l 1.02:1.12
  ti=18 tf=28 e2=1 PlotTI='8 8' FitZVEnergy h${Heavy} h${Heavy} 1.646:1.654
  ti=6  tf=24 ti2=12 e2=1 PlotTI='5 5' FitZVEnergy s s .37:0.42
  ti=9  tf=20 ti2=6  e=1 PlotTI='5 5' PlotTF='26 26' FitZVEnergy l l .17:.21
  #ti=6  tf=20 ti2=6  e2=1 PlotTI='5 5' PlotTF='26 26' FitZVEnergy l l .17:.21
fi

if [ -v DoZV ]; then MakeZV; fi

if [ -v DoPlot ]; then PlotZV; fi

# Original ZV fits from constant
function FitConst()
{
(
  export UnCorr=
  ti=9 tf=14 FitZV h${Heavy} 24 s # Preferred
  ti=5 tf=13 FitZV l 16 s # Preferred ZV_l
  unset ZVFit # End of preferred fit choices
  ti=9 tf=11 FitZV h${Heavy} 20 s
  ti=10 tf=18 FitZV h${Heavy} 28 s
  ti=8 tf=11 FitZV h${Heavy} 20 l # AltZV

  ti=3 tf=9 FitZV l 12 s
  ti=5 tf=15 FitZV l 20 s
  ti=7 tf=17 FitZV l 24 s
  ti=3 tf=9 FitZV l 12 l
  ti=3 tf=13 FitZV l 16 l
  ti=5 tf=15 FitZV l 20 l
  ti=5 tf=17 FitZV l 24 l
)
}

function FitJan24Light()
{
(
  Q=l
  Spec=s
  yrange=0.7:0.73
  e=3 ti='6 8 10' tf='10 12 14' FitZVNew 16 20 24 # Preferred
  unset ZVFit # End of preferred fit choices
  e=3 ti='6 8 10' tf='10 12 14' UnCorr= FitZVNew 16 20 24
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
  Q=h${Heavy}
  Spec=s
  yrange=0.98:1.06
  e=3 ti='8 9' tf='12 15' FitZVNew 20 24
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
    ZVLink=ZVAltZV.txt # Alternate ZV choices are actually the original choices on C1
  done
  )
fi
