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

if [ -v DoFit ]; then
(
  [ -e $ZVFit ] && unset ZVFit #Don't overwrite fit selection if it already exists
  export UnCorr=
  ti=5 tf=7 FitZV h${Heavy} 12 l
  ti=7 tf=9 FitZV h${Heavy} 16 l
  ti=8 tf=12 FitZV h${Heavy} 20 l # AltZV
  ti=8 tf=16 FitZV h${Heavy} 24 l
  ti=8 tf=12 FitZV h${Heavy} 20 s
  ti=8 tf=16 FitZV h${Heavy} 24 s # Preferred
  ti=9 tf=19 FitZV h${Heavy} 28 s
  ti=8 tf=24 FitZV h${Heavy} 32 s

  ti=5 tf=12 FitZV l 16 s # Preferred
  ti=7 tf=18 FitZV l 24 s
  ti=5 tf=13 FitZV l 16 l
  ti=7 tf=19 FitZV l 24 l
)
fi
