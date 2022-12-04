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
MesonFit[s_l,3]=corr_8_28_8_28
MesonFit[s_l,4]=corr_10_28_8_28
MesonFit[s_l,5]=corr_10_28_12_28
MesonFit[s_l,6]=corr_10_28_12_28

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

if [ "${DoAll+x}" = x ]; then
NumExp=2 DeltaT="24 28 32" TI='13 14 14' TF='16 20 24' yrangeR3=0.00068:0.000782 yrangeMEL=0.527:0.537 DoFit 0

NumExp=2 DeltaT="28 32" TI='12 12' TF='20 24' yrangeR3=0.00060:0.00070 yrangeMEL=0.486:0.496 DoFit 1 # Alternate
NumExp=2 DeltaT="24 28 32" TI='13 13 13' TF='16 19 23' yrangeR3=0.00060:0.00070 yrangeMEL=0.486:0.496 DoFit 1 # Preferred

NumExp=2 DeltaT="24 28 32" TI='14 15 15' TF='17 20 24' yrangeR3=0.000185:0.000225 yrangeMEL=0.45:0.46 Gamma=gXYZ DoFit 1

NumExp=2 DeltaT="24 28 32" TI='12 13 13' TF='15 19 23' yrangeR3=0.00056:0.00066 yrangeMEL=0.452:0.462 DoFit 2
NumExp=2 DeltaT="24 28 32" TI='13 13 14' TF='16 20 24' yrangeR3=0.000155:0.000195 yrangeMEL=0.44:0.467 Gamma=gXYZ DoFit 2

NumExp=3 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' yrangeR3=0.000540:0.000625 yrangeMEL=0.42:0.44 DoFit 3
NumExp=2 DeltaT="24 28 32" TI='12 13 14' TF='16 20 24' yrangeR3=0.000130:0.000170 yrangeMEL=0.41:0.44 Gamma=gXYZ DoFit 3

NumExp=2 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' yrangeR3=0.0005:0.00057 yrangeMEL=0.395:0.415 DoFit 4
NumExp=2 DeltaT="24 28" TI='12 13 14' TF='16 20 24' yrangeR3=0.000115:0.00014 yrangeMEL=0.38:0.41 Gamma=gXYZ DoFit 4

NumExp=3 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' yrangeR3=0.00045:0.0006 yrangeMEL=0.36:0.4 DoFit 5
NumExp=2 DeltaT="24 28" TI='12 14 14' TF='16 20 24' yrangeR3=0.00011:0.00014 yrangeMEL=0.34:0.4 Gamma=gXYZ DoFit 5

NumExp=2 DeltaT="24 28 32" TI='12 15 18' TF='15 19 22' yrangeR3=0.00044:0.00052 yrangeMEL=0.32:0.37 DoFit 6
NumExp=2 DeltaT="24 28 32" TI='10 12 16' TF='16 20 24' yrangeR3=0.00008:0.00014 yrangeMEL=0.27:0.37 Gamma=gXYZ DoFit 6
fi
