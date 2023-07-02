#!/usr/bin/env bash

export Ensemble=${Ensemble:-M1}
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

if [ -v DoFit ]; then
  [ -e $ZVFit ] && unset ZVFit #Don't overwrite fit selection if it already exists
  ti=7 tf=9 FitZV h${Heavy} 16 s
  ti=8 tf=12 FitZV h${Heavy} 20 s # Preferred ZV_h
  ti=8 tf=16 FitZV h${Heavy} 24 s
  ti=6 tf=10 FitZV h${Heavy} 16 l
  ti=8 tf=12 FitZV h${Heavy} 20 l
  ti=8 tf=16 FitZV h${Heavy} 24 l

  ti=7 tf=9 FitZV l 16 s
  ti=8 tf=12 FitZV l 20 s
  ti=8 tf=16 FitZV l 24 s # Preferred ZV_l
  ti=6 tf=10 FitZV l 16 l
  ti=8 tf=12 FitZV l 20 l
  ti=8 tf=16 FitZV l 24 l
fi
