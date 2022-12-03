#!/usr/bin/env bash

export Ensemble=F1M
. PlotCommon.sh

#set -x
set -e

############################################################

# Two point fit choices

############################################################

declare -A MesonFit

MesonFit[h385_s,0]=corr_6_29_5_29
MesonFit[s_l,0]=corr_10_26_7_26 # Original with Peter
MesonFit[s_l,1]=corr_8_28_9_28
MesonFit[s_l,2]=corr_8_26_8_28

############################################################

# Perform fit

############################################################

# Parameters:
#   1: pSnk
function DoFit()
{
  local pSnk=$1
  local MesonSnk; GetMesonFile MesonSnk $qSnk $qSpec
  local MesonSrc; GetMesonFile MesonSrc $qSrc $qSpec
  local FitSnk=${MesonFit[$MesonSnk,$pSnk]}
  local FitSrc=${MesonFit[$MesonSrc,0]}
  ( . FitMEL.sh ) #export everything including arrays, but don't upset my environment
}

############################################################

# Now make three-point fits

############################################################

FitWhat="R3"
Ratio=ratioE1ZV1

qSrc=h385
qSnk=l
qSpec=s
Gamma=gT

#YRangeMEL=(0.527:0.537 0.486:0.496 0.452:0.462 0.42:0.44)
#YRangeR3Raw=(0.00068:0.000782 0.00060:0.00070 0.00056:0.00066 0.00052:0.00062)
#YRangeR3=(0.8:1.07 0.75:0.97 0.75:0.92 0.65:0.92)

NumExp=2 DeltaT="24 28 32" TI='13 14 14' TF='16 20 24' yrangeR3=0.00068:0.000782 yrangeMEL=0.527:0.537 DoFit 0

NumExp=2 DeltaT="28 32" TI='12 12' TF='20 24' yrangeR3=0.00060:0.00070 yrangeMEL=0.486:0.496 DoFit 1 # Alternate
NumExp=2 DeltaT="24 28 32" TI='13 13 13' TF='16 19 23' yrangeR3=0.00060:0.00070 yrangeMEL=0.486:0.496 DoFit 1 # Preferred

NumExp=2 DeltaT="24 28 32" TI='14 15 15' TF='17 20 24' yrangeR3=0.000185:0.000225 yrangeMEL=0.45:0.46 Gamma=gXYZ DoFit 1

NumExp=2 DeltaT="24 28 32" TI='12 13 13' TF='15 19 23' yrangeR3=0.00056:0.00066 yrangeMEL=0.452:0.462 DoFit 2
NumExp=2 DeltaT="24 28 32" TI='13 13 14' TF='16 20 24' yrangeR3=0.000155:0.000195 yrangeMEL=0.44:0.467 Gamma=gXYZ DoFit 2
