#!/usr/bin/env bash

#set -x
set -e

#export Ensemble=F1M
export Ratio=ratioE1ZV1

export qSrc=h385
export qSnk=l
export qSpec=s
export pSnk

YRangeMEL=(0.527:0.537 0.486:0.496 0.452:0.462 0.42:0.44)
YRangeR3Raw=(0.00068:0.000782 0.00060:0.00070 0.00056:0.00066 0.00052:0.00062)
YRangeR3=(0.8:1.07 0.75:0.97 0.75:0.92 0.65:0.92)

#export TI='12 13 13 12'
#export TF='15 18 22 26'
export TI='13 15 15'
export TF='16 20 24'
export DeltaT="24 28 32"
export FitWhat="R3"
#export Simul=

# 2pt fit range Source (Heavy) Ds
export Exp2ptSrc=${Exp2ptSrc:-3}
case $Exp2ptSrc in
  2) TISrc=12; TFSrc=29;;
  3) TISrc=6;  TFSrc=29;;
  *) echo "Error Exp2ptSrc=${Exp2ptSrc}"; exit 1;;
esac
export TISrc TFSrc

# 2pt fit range Sink (Light) K
export Exp2ptSnk=${Exp2ptSnk:-3}
case $Exp2ptSnk in
  2) TISnk=7; aTFSnk=(18 18 23 18);;
  3) TISnk=4; aTFSnk=(18 18 23 18);;
  *) echo "Error Exp2ptSnk=${Exp2ptSnk}"; exit 1;;
esac
export TISnk

function PlotOne()
{
  local pSnk=$1
  yrangeMEL=${YRangeMEL[pSnk]} yrangeR3="${YRangeR3Raw[pSnk]}" \
    TFSnk=${aTFSnk[pSnk]} FitMEL.sh
}

#yrange=0.725:0.735 Meson=h385_s p=0 NumExp=3 TI=8 TF=25 LabelTF=35 PlotTDTwoPoint.sh
#yrange='0.185:0.20' Meson=s_l p=0 NumExp=3 TI=4 TF=17 LabelTF=35 PlotTDTwoPoint.sh

PlotOne 0
