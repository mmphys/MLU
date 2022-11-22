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
YRangeR3=(0.8:1.07 0.75:0.97)
YRangeR3Raw=(0.00055:0.0008 *:*)

#export TI='12 13 13 12'
#export TF='15 18 22 26'
export TI='13 15 15'
export TF='16 20 24'
export DeltaT="24 28 32"
export FitWhat="R3"

# Tobi asked for this range
export Exp2pt=${Exp2pt:-3}
#export Simul=
case $Exp2pt in
  2)
  export TISnk=7
  export TFSnk=18
  export TISrc=12
  export TFSrc=29;;
  3)
  export TISnk=4
  export TFSnk=18
  export TISrc=6
  export TFSrc=29;;
esac

#yrange=0.725:0.735 Meson=h385_s p=0 NumExp=3 TI=8 TF=25 LabelTF=35 PlotTDTwoPoint.sh
#yrange='0.185:0.20' Meson=s_l p=0 NumExp=3 TI=4 TF=17 LabelTF=35 PlotTDTwoPoint.sh

pSnk=0; yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} yrangeR3Raw=${YRangeR3Raw[pSnk]} FitMEL.sh
#pSnk=1; yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} yrangeR3Raw=${YRangeR3Raw[pSnk]} FitMEL.sh
