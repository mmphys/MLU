#!/usr/bin/env bash

export Ensemble=${Ensemble:-M3}
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
  ti=10 tf=28 PlotTI='10 10' FitZVEnergy h${Heavy} h${Heavy} 1.229:1.262
  ti=10 tf=24 ti2=18 e2=1 PlotTI='10 11' FitZVEnergy h${Heavy} s 0.815:0.845
  ti=9 tf=24 ti2=14 e2=1 PlotTI='9 9' FitZVEnergy h${Heavy} l 0.787:0.825
  ti=7  tf=20 ti2=13 e2=1 PlotTI='6 6' FitZVEnergy s s .28:.34
  ti=7 tf=18 ti2=10 e2=1 PlotTI='6 5' FitZVEnergy s l 0.22:0.28
  ti=6 tf=18 e2=1 PlotTI='5 5' FitZVEnergy l l .16:.22
fi

if [ -v DoZV ]; then MakeZV; fi

if [ -v DoPlot ]; then PlotZV; fi

if [ -v DoTest ]; then
(
  unset ZVEFit # Don't overwrite selection of fits
  PlotTF='29 29'
)
fi

if [ -v DoFit ]; then
(
  [ -e $ZVFit ] && unset ZVFit #Don't overwrite fit selection if it already exists
  export UnCorr=
  ti=8 tf=12 FitZV h${Heavy} 20 l
  ti=8 tf=14 FitZV h${Heavy} 24 l # Preferred
  ti=8 tf=20 FitZV h${Heavy} 28 l

  ti=4 tf=16 FitZV l 20 s # Preferred
  ti=4 tf=20 FitZV l 24 s
  ti=8 tf=20 FitZV l 28 s
  ti=4 tf=16 FitZV l 20 l
  ti=4 tf=20 FitZV l 24 l
  ti=4 tf=24 FitZV l 28 l
)
fi
