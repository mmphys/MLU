#!/usr/bin/env bash

export Ensemble=${Ensemble:-C2}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
else
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

aMesonFit[h${Heavy}_s,0]=corr_8_24_13_24
aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
aMesonFileMom[h${Heavy}_s,0]=_p2_0

#for FileSeries in ${series-old disp priorPW betterPW priorP betterP}; do
for FileSeries in ${series-disp}; do
  case $FileSeries in
    old) # Versions using different PP+PW fit on each momentum
    aMesonFit[s_l,0]=corr_5_20_7_20
    aMesonFit[s_l,1]=corr_6_22_7_22
    aMesonFit[s_l,2]=corr_5_20_7_20
    aMesonFit[s_l,3]=corr_6_20_7_20
    aMesonFit[s_l,4]=corr_5_18_5_18
    for((i = 0; i < 5; ++i)); do
      aMesonFileOp[s_l,$i]=g5P_g5W
      aMesonFileMom[s_l,$i]=_p2_$i
    done;;

    disp) # Simultaneous fit to PP at all momenta using dispersion relation
    for((i = 0; i < 5; ++i)); do
      aMesonFit[s_l,$i]=corr_5_20_6_22_5_20_6_20_5_18
      aMesonFileOp[s_l,$i]=g5P
      aMesonFileMom[s_l,$i]=
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

  #NumExp=2 DeltaT="20 24 28" TI='8 8 8' TF='15 19 23' yrangeR3=0.0015:0.0019 yrangeMEL=0.77:0.785 DoFit 0
  #NumExp=2 DeltaT="20 24 28" TI='9 9 9' TF='13 16 20' yrangeR3=0.0014:0.00165 yrangeMEL=0.68:0.71 DoFit 1

done
fi
#
  NumExp=2 DeltaT="20 24 28" TI='8 8 8' TF='14 18 22' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 DoFit 2
