#!/usr/bin/env bash

export Ensemble=${Ensemble:-C2}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit 2
fi
. FitMEL.sh

#set -x
set -e

############################################################

# Two point fit choices

############################################################

# y axis ranges for 2pt fit plots
declare -A ayRange
ayRange[h6413_s,0]='1.095:1.135'
ayRange[s_l,0]='0.315:0.365'
ayRange[s_l,1]='0.4:0.45'
ayRange[s_l,2]='0.45:0.56'
ayRange[s_l,3]='0.51:0.64'
ayRange[s_l,4]='0.57:0.75'

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
ayrangeR3[gXYZ,4]=0.00020:0.00040
ayrangeMEL[gXYZ,1]=0.6:0.76
ayrangeMEL[gXYZ,2]=0.48:0.66
ayrangeMEL[gXYZ,3]=0.46:0.62
ayrangeMEL[gXYZ,4]=0.42:0.58

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
aMesonFileMom[h${Heavy}_s,0]=_p2_0

############################################################

# Choose two step fits - i.e. fit R3 using 2pt fits as input

############################################################

function ChooseTwoPtFits()
{
  unset FitOptions
  local DsTIP=8
  local DsTFP=24
  local DsTIW=13
  local DsTFW=24
  aMesonFit[h${Heavy}_s,0]=corr_${DsTIP}_${DsTFP}_${DsTIW}_${DsTFW}

  case "$1" in
    old) # Versions using different PP+PW fit on each momentum
    local aKaonTIP=( 5  6  5  6  5)
    #local aKaonTFP=(20 21 20 19 18)
    local aKaonTFP=(20 22 19 20 17) # 16 May 2023
    local aKaonTIW=( 7  7  7  7  5)
    local aKaonTFW=(20 22 20 20 18)
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
      # aMesonFit[s_l,$i]=corr_5_20_6_21_5_20_6_19_5_18 # thinned 1:3:2
      aMesonFit[s_l,$i]=corr_5_20_6_22_5_19_6_20_5_17 # thinned 2
      aMesonFileOp[s_l,$i]=g5P
      aMesonFileMom[s_l,$i]=
    done;;

    dispind) # Simultaneous fit to PP at all momenta using dispersion relation
    FitOptions='--nopolap g5P'
    aMesonFit[h${Heavy}_s,0]=$FileSeries.corr_${DsTIP}_${DsTFP}_${DsTIW}_${DsTFW}
    for((i = 0; i < 5; ++i)); do
      aMesonFit[s_l,$i]=$FileSeries.corr_5_20_6_22_5_19_6_20_5_17 # thinned 2
      aMesonFileOp[s_l,$i]=g5P
      aMesonFileMom[s_l,$i]=
    done;;

    *) echo "Two-point fits $1 unrecognised"; exit 1;;
  esac
}

############################################################

# These were my original ratio fit ranges

############################################################

function RatioFitsBase()
{
    Gamma=gT
    #Thinning="'' 2 3"
    NumExp=3 DeltaT="16 20 24" TI='9 10 10' TF='10 13 17' FitTwoStage 0

    NumExp=3 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 1
    #NumExp=3 DeltaT="20 24" TI='10 10' TF='14 18' FitTwoStage 1 #Bad
    #NumExp=3 DeltaT="20 24" TI='9 9' TF='14 18' FitTwoStage 1 #Bad
    #NumExp=3 DeltaT="20 24 28" TI='9 9 9' TF='14 18 22' FitTwoStage 1 #Bad
    #NumExp=3 DeltaT="16 20 24" TI='9 10 10' TF='10 13 17' FitTwoStage 1 #Bad
    #NumExp=3 DeltaT="16 20 24" TI='9 9 9' TF='10 13 16' FitTwoStage 1 # Bad
    #NumExp=3 DeltaT="16 20 24" TI='9 9 9' TF='10 14 18' FitTwoStage 1 #Bad
    #NumExp=3 DeltaT="16 20 24" TI='8 8 8' TF='10 14 18' FitTwoStage 1 #Bad

    NumExp=3 DeltaT="16 20" TI='8 8' TF='10 13' FitTwoStage 2
    #NumExp=3 DeltaT="16 20 24" TI='9 10 10' TF='10 13 17' FitTwoStage 2 #Bad
    #NumExp=3 DeltaT="16 20" TI='9 10' TF='10 13' FitTwoStage 2 # Good
    #NumExp=3 DeltaT="16 20" TI='8 10' TF='10 13' FitTwoStage 2 # Good

    NumExp=3 DeltaT="16 20" TI='8 8' TF='10 14' FitTwoStage 3
    #NumExp=3 DeltaT="16 20 24" TI='8 8 8' TF='10 13 17' FitTwoStage 3 # Very Bad
    #NumExp=3 DeltaT="16 20" TI='8 8' TF='10 13' FitTwoStage 3 # OK

    NumExp=3 DeltaT="16 20" TI='9 9' TF='10 14' FitTwoStage 4 #Best
    #NumExp=3 DeltaT="16 20" TI='8 9' TF='10 14' FitTwoStage 4 #Almost as good
    #NumExp=3 DeltaT="16 20" TI='8 8' TF='10 14' FitTwoStage 4 #Good
    #NumExp=3 DeltaT="16 20" TI='7 8' TF='10 14' FitTwoStage 4 #Bad
    #NumExp=3 DeltaT="16 20" TI='7 8' TF='11 15' FitTwoStage 4 #Bad
    #NumExp=3 DeltaT="12 16 20" TI='6 8 8' TF='7 10 14' FitTwoStage 4 #Bad

    Gamma=gXYZ
    NumExp=3 DeltaT="16 20 24" TI='9 9 9' TF='11 15 19' FitTwoStage 1 #Good
    #NumExp=3 DeltaT="16 20 24" TI='8 8 8' TF='10 14 18' FitTwoStage 1 #Good
    #NumExp=3 DeltaT="16 20 24" TI='9 9 12' TF='11 15 19' FitTwoStage 1 #Good
    #NumExp=3 DeltaT="16 20 24" TI='9 9 10' TF='11 15 19' FitTwoStage 1 #Good

    NumExp=3 DeltaT="16 20 24" TI='8 9 9' TF='10 14 18' FitTwoStage 2 #Good
    #NumExp=3 DeltaT="16 20"    TI='8 9'   TF='10 14'    FitTwoStage 2 #Good
    #NumExp=3 DeltaT="16 20 24" TI='9 9 9' TF='11 15 19' FitTwoStage 2 #Good
    #NumExp=3 DeltaT="16 20 24" TI='8 9 9' TF='10 15 19' FitTwoStage 2 #Good

    NumExp=3 DeltaT="16 20 24" TI='8 9 9' TF='10 14 18' FitTwoStage 3 #Good
    #NumExp=3 DeltaT="16 20"    TI='8 9'   TF='10 14'    FitTwoStage 3 #Good
    #NumExp=3 DeltaT="16 20 24" TI='8 9 12' TF='10 14 18' FitTwoStage 3 #Good
    #NumExp=3 DeltaT="16 20 24" TI='8 9 14' TF='10 14 18' FitTwoStage 3 #Good

    NumExp=3 DeltaT="16 20" TI='7 9' TF='11 14' FitTwoStage 4 #Good
    #NumExp=3 DeltaT="12 16 20" TI='6 7 9' TF='7 11 14' FitTwoStage 4 #Interesting - believably higher
    #NumExp=3 DeltaT="16 20 24" TI='8 9 9' TF='10 14 18' FitTwoStage 4 #Bad
    #NumExp=3 DeltaT="16 20" TI='8 9' TF='10 14' FitTwoStage 4 #So-So
    #NumExp=3 DeltaT="16 20" TI='8 9' TF='11 14' FitTwoStage 4 #OK
}

############################################################

# Simultaneously fit two-pt functions and ratios

############################################################

function RatioFitsSimul()
{
  Gamma=gT
  DeltaT=(16 20 24)
  TI=(8 8 8)
  TF=(10 14 20)
  Thin=('' '' t3)
  NumExp=2
  for (( n=$MaxPSq; n<=$MaxPSq; ++n )); do
    if ((n >= 2)); then DeltaT=(16 20); fi
    DoSimulFit $n
    if ((n)); then
      Gamma=gXYZ DoSimulFit $n
      IncludeSpatial= DoSimulFit $n
    fi
  done
}

############################################################

# Ratio fit ranges being tested

############################################################

function RatioFitsTest()
{
  for (( i=0; i<2; ++i )); do
  Gamma=gT
  NumExp=2
  Ratio=ratio; Renorm=; NotRaw=
  #DeltaT="20 24 28" TI='8 8 8' TF='15 19 23' yrangeR3=0.0015:0.0019 yrangeMEL=0.77:0.785 FitTwoStage 0
  #DeltaT="20 24 28" TI='9 9 9' TF='13 16 20' yrangeR3=0.0014:0.00165 yrangeMEL=0.68:0.71 FitTwoStage 1
  #NumExp=2 DeltaT="20 24" TI='8 8 9' TF='14 17 21' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="2 3 4" FitTwoStage 2 # As at 18 Apr 2023

  #NumExp=2 DeltaT="20 24" TI='8 8 9' TF='14 17 21' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="3 3 4" Thinning_Save="2:1:1:2:2:1:1 3 4"  FitTwoStage 2

  #NumExp=3 DeltaT="16 20 24" TI='8 8 8' TF='10 15 17' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="'' '' 3" FitTwoStage 2 # THIS IS THE BEST FIT p=0.05

  # Kick out DT=16 for Tobi - hmmm ... this looks even better!
  #NumExp=2 DeltaT="20 24" TI='8 8' TF='15 17' yrangeR3=0.0012:0.0016 yrangeMEL=0.59:0.65 Thinning="'' 3" FitTwoStage 2 # THIS IS THE BEST FIT p=0.05

  # With Peter 19 Apr
  #NumExp=2 DeltaT="16 20 24" TI='8 8 8' TF='10 14 18' yrangeR3=0.00115:0.0016 yrangeMEL=0.59:0.65 Thinning="" FitTwoStage 2 # THIS IS THE BEST FIT p=0.05
  UnCorr=
  done
}

############################################################

# Now make three-point fits

############################################################

FitWhat="${FitWhat-R3}"
Ratio=ratioE1ZV1

qSrc=h$Heavy
qSnk=l
qSpec=s

for FileSeries in ${series-disp dispind renorm}; do
(
  echo "Performing $FitWhat fits to $Ratio for $FileSeries"
  case $FileSeries in
    old)
      ChooseTwoPtFits $FileSeries
      RatioFitsBase;;

    disp | dispind)
      ChooseTwoPtFits $FileSeries
      RatioFitsBase
      [ -v UnCorr ] || UnCorr= RatioFitsBase;;

    renorm)
      ChooseTwoPtFits disp
      Ratio=ratio Renorm= NotRaw= RatioFitsBase;;

    simul)
      RatioFitsSimul;;

    "test")
      ChooseTwoPtFits disp
      RatioFitsTest;;

    *) echo "FileSeries $FileSeries unrecognised"; exit 1;;
  esac
)
done
