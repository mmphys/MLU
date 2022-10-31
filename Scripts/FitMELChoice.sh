#!/usr/bin/env bash

#set -x
set -e

export qSrc=h385
export qSnk=l
export qSpec=s
export pSnk

YRangeMEL=(0.524:0.539 0.48:0.496)
YRangeR3=(0.75:1.05 0.75:0.97)

#pSnk=0 TI='11 12 14 12' TF='15 19 21 26' yrangeMEL=0.52:0.537 yrangeR3=0.75:1.05 FitMEL.sh
#pSnk=1 TI='10 12 14 10' TF='13 16 20 26' yrangeMEL=0.48:0.496 yrangeR3=0.7:0.96 FitMEL.sh

# Tobi asked for this range
#export Exp2pt=3
#export TISnk=9
#export TFSnk=22
#export TISrc=8
#export TFSrc=23
#export Simul=1

export Exp2pt=2
export TISnk=11
export TFSnk=22
export TISrc=10
export TFSrc=23
export Simul=1

export TI='11 11 13 13'
export TF='14 18 22 26'
pSnk=0; yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} FitMEL.sh
#pSnk=1; yrangeMEL=${YRangeMEL[pSnk]} yrangeR3=${YRangeR3[pSnk]} FitMEL.sh
