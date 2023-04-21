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

declare -A ayrangeR3
declare -A ayrangeMEL
ayrangeR3[gT,0]=0.0016:0.00185
ayrangeR3[gT,1]=0.0014:0.00165
ayrangeR3[gT,2]=0.0011:0.0016
ayrangeR3[gT,3]=0.0010:0.0015
ayrangeR3[gT,4]=0.0008:0.0013
ayrangeMEL[gT,0]=0.77:0.785
ayrangeMEL[gT,1]=0.68:0.71
ayrangeMEL[gT,2]=0.59:0.65
ayrangeMEL[gT,3]=0.52:0.58
ayrangeMEL[gT,4]=0.47:0.55
ayrangeR3[gXYZ,1]=0.00040:0.0006
ayrangeR3[gXYZ,2]=0.00030:0.0005
ayrangeR3[gXYZ,3]=0.00025:0.00045
ayrangeR3[gXYZ,4]=0.00015:0.00035
ayrangeMEL[gXYZ,1]='*:*'
ayrangeMEL[gXYZ,2]='*:*'
ayrangeMEL[gXYZ,3]='*:*'
ayrangeMEL[gXYZ,4]='*:*'

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

aMesonFit[h${Heavy}_s,0]=corr_8_24_13_24
aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
aMesonFileMom[h${Heavy}_s,0]=_p2_0

aKaonTIP=(5 6 5 6 5)
aKaonTFP=(20 22 20 20 18)
aKaonTIW=(7 7 7 7 5)
aKaonTFW=(20 22 20 20 18)

#for FileSeries in ${series-old disp priorPW betterPW priorP betterP}; do
for FileSeries in ${series-disp}; do
  case $FileSeries in
    old) # Versions using different PP+PW fit on each momentum
    for((i = 0; i < ${#aKaonTIP[@]}; ++i)); do
      aMesonFit[s_l,$i]=corr_${aKaonTIP[i]}_${aKaonTFP[i]}_${aKaonTIW[i]}_${aKaonTFW[i]}
    done
    for((i = 0; i < 5; ++i)); do
      aMesonFileOp[s_l,$i]=g5P_g5W
      aMesonFileMom[s_l,$i]=_p2_$i
    done;;

    disp) # Simultaneous fit to PP at all momenta using dispersion relation
    for((i = 0; i < 5; ++i)); do
      # aMesonFit[s_l,$i]=corr_5_20_6_22_5_20_6_20_5_18 # unthinned
      aMesonFit[s_l,$i]=corr_5_20_6_21_5_20_6_19_5_18 # thinned 1:3:2
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

  (
    UnCorr=
    DeltaT="16 20 24"
    TI='8 8 8'
    TF='10 14 18'
    NumExp=2
    for (( n=0; n<5; ++n )); do
      DoFit $n
      if ((n)); then Gamma=gXYZ DoFit $n; fi
    done
  )
done
fi
#
  #NumExp=2 DeltaT="20 24" TI='8 8 9' TF='14 17 21' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="2 3 4" DoFit 2 # As at 18 Apr 2023

  #NumExp=2 DeltaT="20 24" TI='8 8 9' TF='14 17 21' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="3 3 4" Thinning_Save="2:1:1:2:2:1:1 3 4"  DoFit 2

  #NumExp=3 DeltaT="16 20 24" TI='8 8 8' TF='10 15 17' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="'' '' 3" DoFit 2 # THIS IS THE BEST FIT p=0.05

  # Kick out DT=16 for Tobi - hmmm ... this looks even better!
  #NumExp=2 DeltaT="20 24" TI='8 8' TF='15 17' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="'' 3" DoFit 2 # THIS IS THE BEST FIT p=0.05

  # With Peter 19 Apr
  #NumExp=2 DeltaT="16 20 24" TI='8 8 8' TF='10 14 18' yrangeR3=0.00115:0.0016 yrangeMEL=0.59:0.65 Thinning="" DoFit 2 # THIS IS THE BEST FIT p=0.05
