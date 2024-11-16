#!/usr/bin/env bash

# Ensemble parameters
Ensemble=F1M
L=48
T=96
EnsembleaInv=2.708 #GeV

# Maximum psquared
MaxPSq=6

EnsembleDeltaT=(16 20 24 28 32)
Heavy=385

export MLUSeed=3285139232

[ -e "$HOME/.MLUConfig.sh" ] && . "$HOME/.MLUConfig.sh"
[ -z "$PlotData" ] && PlotData=/Volumes/QCD/tursa/semilep/data/$Ensemble/analyse
