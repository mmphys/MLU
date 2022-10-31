#!/usr/bin/env bash

#set -x
set -e

export Ensemble=F1M
export Ratio=ratio.manual

export qSrc=h385
export qSnk=l
export qSpec=s
export pSnk

YRangeMEL=(0.524:0.539 0.48:0.496)
YRangeR3=(0.75:1.05 0.75:0.97)

export TI='11 11 11 11'
export TF='14 18 22 26'

# Tobi asked for this range
export Exp2pt=2
#export Simul=
case $Exp2pt in
  2)
  export TISnk=9
  export TFSnk=22
  export TISrc=12
  export TFSrc=25;;
  3)
  export TISnk=9
  export TFSnk=22
  export TISrc=8
  export TFSrc=23;;
esac

for (( pSnk=0; pSnk <=6; ++pSnk ))
do
  yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} FitMEL.sh
done
