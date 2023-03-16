#!/usr/bin/env bash

export Ensemble=C1
. PlotCommon.sh

#set -x
set -e

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
  local FitSnk=${aMesonFit[$MesonSnk,$pSnk]}
  local FitSrc=${aMesonFit[$MesonSrc,0]}
  local FileOpSnk=${aMesonFileOp[$MesonSnk,$pSnk]}
  local FileOpSrc=${aMesonFileOp[$MesonSrc,0]}
  local FileMomSnk=${aMesonFileMom[$MesonSnk,$pSnk]}
  local FileMomSrc=${aMesonFileMom[$MesonSrc,0]}
  ( . FitMEL.sh ) #export everything including arrays, but don't upset my environment
}

############################################################

# Two point fit choices

############################################################

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

aMesonFit[h${Heavy}_s,0]=corr_6_27_14_27
aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
aMesonFileMom[h${Heavy}_s,0]=_p2_0

for FileSeries in ${series-old disp priorPW priorP}; do
case $FileSeries in
  old) # Versions using different PP+PW fit on each momentum
  aMesonFit[s_l,0]=corr_6_23_7_23
  aMesonFit[s_l,1]=corr_6_23_7_23
  aMesonFit[s_l,2]=corr_5_20_5_20
  aMesonFit[s_l,3]=corr_5_20_5_20
  aMesonFit[s_l,4]=corr_5_18_5_18
  for((i = 0; i < 5; ++i)); do
    aMesonFileOp[s_l,$i]=g5P_g5W
    aMesonFileMom[s_l,$i]=_p2_$i
  done;;

  disp) # Simultaneous fit to PP at all momenta using dispersion relation
  for((i = 0; i < 5; ++i)); do
    aMesonFit[s_l,$i]=corr_6_23_6_23_5_20_5_20_5_18
    aMesonFileOp[s_l,$i]=g5P
  done;;

  priorPW) # Simultaneous fit to PP & PW at all momenta using dispersion relation
  aMesonFit[s_l,0]=corr_6_23_7_23
  aMesonFit[s_l,1]=corr_6_23_7_23
  aMesonFit[s_l,2]=corr_5_20_5_20
  aMesonFit[s_l,3]=corr_5_20_5_20
  aMesonFit[s_l,4]=corr_5_18_5_18
  aMesonFileOp[s_l,0]=g5P_g5W
  aMesonFileMom[s_l,0]=_p2_0
  for((i = 1; i < 5; ++i)); do
    aMesonFileOp[s_l,$i]=g5P_g5W
    aMesonFileMom[s_l,$i]=_p2_$i.${FileSeries}_6_23_7_23
  done;;

  priorP) # Simultaneous fit to PP at all momenta using dispersion relation
  aMesonFit[s_l,0]=corr_6_23_7_23
  aMesonFit[s_l,1]=corr_6_23
  aMesonFit[s_l,2]=corr_5_20
  aMesonFit[s_l,3]=corr_5_20
  aMesonFit[s_l,4]=corr_5_18
  aMesonFileOp[s_l,0]=g5P_g5W
  aMesonFileMom[s_l,0]=_p2_0
  for((i = 1; i < 5; ++i)); do
    aMesonFileOp[s_l,$i]=g5P_g5W
    aMesonFileMom[s_l,$i]=_p2_$i.${FileSeries}_6_23_7_23
  done;;

  *) echo "FileSeries $FileSeries unrecognised"; exit 1;;
esac

############################################################

# Now make three-point fits

############################################################

FitWhat="R3"
Ratio=ratioE1ZV1

qSrc=h$Heavy
qSnk=l
qSpec=s
Gamma=gT

echo "Performing $FitWhat fits to $Ratio for $FileSeries"

  NumExp=2 DeltaT="20 24 28" TI='9 11 13' TF='14 18 18' yrangeR3=0.0015:0.0019 yrangeMEL=0.76:0.8 DoFit 0 # Preferred
  NumExp=2 DeltaT="20 24 28" TI='11 14 13' TF='12 15 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.68:0.71 DoFit 1
  NumExp=2 DeltaT="20 24 28" TI='8 12 16' TF='12 16 19' yrangeR3=0.0013:0.00160 yrangeMEL=0.58:0.66 DoFit 2
  NumExp=2 DeltaT="20 24 28" TI='8 12 16' TF='14 16 19' yrangeR3=0.0010:0.00160 yrangeMEL=0.50:0.64 DoFit 3
  NumExp=2 DeltaT="20 24" TI='8 12' TF='14 16' yrangeR3=0.0010:0.00150 yrangeMEL=0.44:0.62 DoFit 4

  # Gamma spatial
  NumExp=2 DeltaT="20 24 28" TI='12 15 13' TF='14 17 17' yrangeR3=0.0005:0.00065 yrangeMEL=0.66:0.74 Gamma=gXYZ DoFit 1
  NumExp=2 DeltaT="16 20 24" TI='8 9 8' TF='11 15 19' Gamma=gXYZ yrangeR3=0.00038:0.00055 yrangeMEL=0.56:0.70 DoFit 2 # My preferred fit
  NumExp=2 DeltaT="16 20 24" TI='8 9 8' TF='11 15 19' Gamma=gXYZ yrangeR3=0.00032:0.00046 yrangeMEL=0.48:0.65 DoFit 3 # My preferred fit
  NumExp=2 DeltaT="16 20 24" TI='7 9 8 13' TF='11 15 19 16' Gamma=gXYZ yrangeR3=0.00026:0.00046 yrangeMEL=0.44:0.65 DoFit 4

if [ "${DoAlt+x}" = x ]; then
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 11 13' TF='11 14 18 18' yrangeR3=0.0015:0.0019 yrangeMEL=0.76:0.8 DoFit 0 # Very much the same as prior fit
  NumExp=1 DeltaT="20 24 28" TI='11 14 13' TF='12 15 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705 DoFit 1 # I prefer this one
  NumExp=2 DeltaT="16 20 24 28" TI='9 11 13 13' TF='10 13 16 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705 DoFit 1
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 8 13' TF='11 15 19 16' Gamma=gXYZ yrangeR3=0.00038:0.00055 yrangeMEL=0.56:0.70 DoFit 2 # Just to see the DT=28 data
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 8 13' TF='11 15 19 16' Gamma=gXYZ yrangeR3=0.00032:0.00046 yrangeMEL=0.48:0.65 DoFit 3 # Just to see the DT=28 data
fi
done

##NumExp=2 DeltaT="20 24 28" TI='10 10 13' TF='13 18 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705 DoFit 1 # I prefer this one
#NumExp=2 DeltaT="20 24 28" TI='10 14 13' TF='13 15 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705 DoFit 1 # I prefer this one
# This works
#NumExp=2 DeltaT="16 20" TI='6 8 12' TF='9 14 16' yrangeR3=0.0008:0.00150 yrangeMEL=0.44:0.62 DoFit 4
# This works - end
