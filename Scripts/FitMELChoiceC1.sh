#!/usr/bin/env bash

#set -x
set -e

export Ensemble=C1

export qSrc=h6413
export qSnk=l
export qSpec=s
export pSnk

YRangeMEL=(*:*)
YRangeR3=(*:*)

export TI='7  7  7  7'
export TF='16 20 24 28'

# Tobi asked for this range
export Exp2pt=2
#export Simul=
case $Exp2pt in
  2)
  export TISnk=6
  export TFSnk=14
  export TISrc=8
  export TFSrc=16;;
  3)
  export TISnk=6
  export TFSnk=14
  export TISrc=5
  export TFSrc=15;;
esac

pSnk=0; yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} FitMEL.sh
