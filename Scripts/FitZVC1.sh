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

if [ -v DoFit ]; then
(
  [ -e $ZVFit ] && unset ZVFit #Don't overwrite fit selection if it already exists
  export UnCorr=
  ti=9 tf=11 FitZV h${Heavy} 20 s
  ti=9 tf=14 FitZV h${Heavy} 24 s # Preferred
  ti=10 tf=18 FitZV h${Heavy} 28 s
  ti=8 tf=11 FitZV h${Heavy} 20 l # AltZV

  ti=3 tf=9 FitZV l 12 s
  ti=5 tf=13 FitZV l 16 s # Preferred ZV_l
  ti=5 tf=15 FitZV l 20 s
  ti=7 tf=17 FitZV l 24 s
  ti=3 tf=9 FitZV l 12 l
  ti=3 tf=13 FitZV l 16 l
  ti=5 tf=15 FitZV l 20 l
  ti=5 tf=17 FitZV l 24 l
)
fi
