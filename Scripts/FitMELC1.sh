#!/usr/bin/env bash

export Ensemble=C1
. PlotCommon.sh

#set -x
set -e

############################################################

# Two point fit choices

############################################################

declare -A MesonFit

#MesonFit[h${Heavy}_s,0]=corr_6_29_5_29 # Preferred fit
#MesonFit[h${Heavy}_s,0]=corr_6_29_13_27
MesonFit[h${Heavy}_s,0]=corr_6_29_17_27

## TEMP - was using this Thu 19 Jan & Mon 23 Jan
#MesonFit[h${Heavy}_s,0]=corr_6_27_5_27
## TEMP END

MesonFit[s_l,0]=corr_6_22_5_23 # Preferred fit
#MesonFit[s_l,0]=corr_6_22_7_23
#MesonFit[s_l,0]=corr_6_22_12_23

MesonFit[s_l,1]=corr_6_22_5_23 # Preferred fit
#MesonFit[s_l,1]=corr_6_14_5_23
#MesonFit[s_l,1]=corr_6_14_12_23
#MesonFit[s_l,1]=corr_6_22_12_23

MesonFit[s_l,2]=corr_6_22_5_23 # 20221211
MesonFit[s_l,3]=corr_6_22_5_23 # 20221211
MesonFit[s_l,4]=corr_6_22_5_23 # 20221211

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

qSrc=h$Heavy
qSnk=l
qSpec=s
Gamma=gT

if [ "${DoAll+x}" = x ]; then
  NumExp=2 DeltaT="20 24 28" TI='9 11 13' TF='14 18 18' yrangeR3=0.0015:0.0019 yrangeMEL=0.76:0.8 DoFit 0 # Preferred
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 11 13' TF='11 14 18 18' yrangeR3=0.0015:0.0019 yrangeMEL=0.76:0.8 DoFit 0 # Very much the same as prior fit
  # Various choices for n^2=1
  yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705
  NumExp=1 DeltaT="20 24 28" TI='11 14 13' TF='12 15 17' DoFit 1 # I prefer this one
  NumExp=2 DeltaT="20 24 28" TI='11 14 13' TF='12 15 17' DoFit 1
  NumExp=2 DeltaT="16 20 24 28" TI='9 11 13 13' TF='10 13 16 17' DoFit 1
  # n^2 >= 2
  NumExp=2 DeltaT="20 24 28" TI='8 12 16' TF='12 16 19' yrangeR3=0.0013:0.00160 yrangeMEL=0.58:0.66 DoFit 2
  NumExp=2 DeltaT="20 24 28" TI='8 12 16' TF='14 16 19' yrangeR3=0.0010:0.00160 yrangeMEL=0.50:0.64 DoFit 3
  NumExp=2 DeltaT="20 24" TI='8 12 16' TF='14 16 19' yrangeR3=0.0010:0.00150 yrangeMEL=0.44:0.62 DoFit 4

  # Gamma spatial
  NumExp=2 DeltaT="20 24 28" TI='12 15 13' TF='14 17 17' yrangeR3=0.0005:0.00065 yrangeMEL=0.66:0.74 Gamma=gXYZ DoFit 1
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 8 13' TF='11 15 19 16' Gamma=gXYZ yrangeR3=0.00038:0.00055 yrangeMEL=0.56:0.70 DoFit 2 # My preferred fit
  NumExp=2 DeltaT="16 20 24" TI='8 9 8' TF='11 15 19' Gamma=gXYZ yrangeR3=0.00038:0.00055 yrangeMEL=0.56:0.70 DoFit 2 # Just to see the DT=28 data
  NumExp=2 DeltaT="16 20 24" TI='8 9 8' TF='11 15 19' Gamma=gXYZ yrangeR3=0.00032:0.00046 yrangeMEL=0.48:0.65 DoFit 3 # My preferred fit
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 8 13' TF='11 15 19 16' Gamma=gXYZ yrangeR3=0.00032:0.00046 yrangeMEL=0.48:0.65 DoFit 3 # Just to see the DT=28 data
  NumExp=2 DeltaT="16 20 24" TI='7 9 8 13' TF='11 15 19 16' Gamma=gXYZ yrangeR3=0.00026:0.00046 yrangeMEL=0.44:0.65 DoFit 4

fi
