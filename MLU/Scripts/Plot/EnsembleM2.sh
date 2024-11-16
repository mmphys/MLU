#!/usr/bin/env bash

# Ensemble parameters
Ensemble=M2
L=32
T=64
EnsembleaInv=2.3833 #GeV

# Maximum psquared
MaxPSq=4

EnsembleDeltaT=(16 20 24 28 32)
Heavy=447

export MLUSeed=2475361470

[ -e "$HOME/.MLUConfig.sh" ] && . "$HOME/.MLUConfig.sh"
[ -z "$PlotData" ] && PlotData=/Volumes/QCD/tursa/semilep/data/$Ensemble/analyse
