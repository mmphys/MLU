#!/usr/bin/env bash

# Ensemble parameters
Ensemble=C1
L=24
T=64
EnsembleaInv=1.7848 #GeV

# Maximum psquared
MaxPSq=4

EnsembleDeltaT=(16 20 24 28 32 12)
Heavy=6413

export MLUSeed=2263212701

[ -e "$HOME/.MLUConfig.sh" ] && . "$HOME/.MLUConfig.sh"
[ -z "$PlotData" ] && PlotData=/Volumes/QCD/tursa/semilep/data/$Ensemble/analyse
