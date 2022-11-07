#!/usr/bin/env bash

#set -x
set -e

export Ensemble=M1

export qSrc=h447
export qSnk=l
export qSpec=s
export pSnk

YRangeMEL=(*:*)
YRangeR3=(*:*)

export TI='10 10 10 10'
export TF='15 19 23 27'

# Tobi asked for this range
export Exp2pt=2
#export Simul=
case $Exp2pt in
  2)
  export TISnk=6
  export TFSnk=23
  export TISrc=9
  export TFSrc=25;;
  2-propto-F1M)
  export TISnk=8
  export TFSnk=19
  export TISrc=11
  export TFSrc=22;;
  3)
  export TISnk=8
  export TFSnk=19
  export TISrc=7
  export TFSrc=20;;
esac

pSnk=0; yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} FitMEL.sh
