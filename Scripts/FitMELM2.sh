#!/usr/bin/env bash

export Ensemble=${Ensemble:-M2}
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
ayRange[h${Heavy}_s,0]='0.815:0.845'
ayRange[s_l,0]='0.22:0.28'
ayRange[s_l,1]='0.28:0.36'
ayRange[s_l,2]='0.32:0.43'
ayRange[s_l,3]='0.34:0.55'
ayRange[s_l,4]='0.37:0.62'

declare -A ayrangeR3
declare -A ayrangeMEL
ayrangeR3[gT,0]=0.0014:0.0019
ayrangeMEL[gT,0]=0.58:0.61
ayrangeR3[gT,1]=0.0012:0.0018
ayrangeMEL[gT,1]=0.51:0.54
ayrangeR3[gT,2]=0.0010:0.0016
ayrangeMEL[gT,2]=0.45:0.50
ayrangeR3[gT,3]=0.0009:0.0015
ayrangeMEL[gT,3]=0.34:0.45
ayrangeR3[gT,4]=0.0008:0.0016
ayrangeMEL[gT,4]=0.30:0.46
ayrangeR3[gXYZ,1]=0.00040:0.00065
ayrangeMEL[gXYZ,1]=0.46:0.56
ayrangeR3[gXYZ,2]=0.00030:0.00055
ayrangeMEL[gXYZ,2]=0.42:0.50
ayrangeR3[gXYZ,3]=0.00025:0.00055
ayrangeMEL[gXYZ,3]=0.36:0.48
ayrangeR3[gXYZ,4]=0.00020:0.00050
ayrangeMEL[gXYZ,4]=0.30:0.48

# y-ranges for renormalised R3 ratios
declare -A ayrangeR3R
ayrangeR3R[gT,0]=0.86:1.1
ayrangeR3R[gT,1]=0.70:0.94
ayrangeR3R[gT,2]=0.62:0.86
ayrangeR3R[gT,3]=0.58:0.82
ayrangeR3R[gT,4]=0.56:0.8
ayrangeR3R[gXYZ,1]=0.25:0.35
ayrangeR3R[gXYZ,2]=0.20:0.30
ayrangeR3R[gXYZ,3]=0.15:0.25
ayrangeR3R[gXYZ,4]=0.12:0.22

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

############################################################

# Choose heavy two point fits

############################################################

function ChooseHeavyFits()
{
  local aDsTIP=10
  local aDsTFP=26
  local aDsTIW=19
  local aDsTFW=26

if ! [ -v DisableDsWall ]; then
  # Ds point-point and point-wall (though they have some tension)
  aMesonFit[h${Heavy}_s,0]=corr_${aDsTIP}_${aDsTFP}_${aDsTIW}_${aDsTFW}
  aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
else
  # Ds point-point only. Is this any better?
  aMesonFit[h${Heavy}_s,0]=corr_${aDsTIP}_${aDsTFP}
  aMesonFileOp[h${Heavy}_s,0]=g5P
fi
aMesonFileMom[h${Heavy}_s,0]=_p2_0
}

############################################################

# Choose two step fits - i.e. fit R3 using 2pt fits as input

############################################################

function ChooseTwoPtFits()
{
  case "$1" in
    old) # Versions using different PP+PW fit on each momentum
    local aKaonTIP=( 6  6  7  7  6)
    local aKaonTFP=(20 20 20 19 16)
    local aKaonTIW=( 8  7  6  6 13)
    local aKaonTFW=(20 20 20 19 16)
    for((i = 0; i < ${#aKaonTIP[@]}; ++i)); do
      aMesonFit[s_l,$i]=corr_${aKaonTIP[i]}_${aKaonTFP[i]}_${aKaonTIW[i]}_${aKaonTFW[i]}
    done
    for((i = 0; i < 5; ++i)); do
      aMesonFileOp[s_l,$i]=g5P_g5W
      aMesonFileMom[s_l,$i]=_p2_$i
    done;;

    disp) # Simultaneous fit to PP at all momenta using dispersion relation
    for((i = 0; i < 5; ++i)); do
      aMesonFit[s_l,$i]=corr_6_20_6_19_7_20_7_18_6_15 # thinned 1:3:2 # Preferred
      #aMesonFit[s_l,$i]=corr_6_20_6_20_7_19_7_19_6_16 # thinned 2
      #aMesonFit[s_l,$i]=corr_6_20_6_20_7_20_7_19_6_16 # unthinned
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
  echo "M2 performing Base ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  # For now, these are borrowed directly from M1
  Gamma=gT
  NumExp=${NumExp:-2}
  DeltaT="24 28 32" TI='11 11 11' TF='18 22 26' FitTwoStage 0
  Alt= DeltaT="28 32" TI='11 11' TF='22 26' FitTwoStage 0
  Alt= DeltaT="20 24 28 32" TI='10 11 11 11' TF='14 18 22 26' FitTwoStage 0
  Alt= DeltaT="20 24 28   " TI='10 11 11'    TF='14 18 22' FitTwoStage 0
  Alt= DeltaT="20 24      " TI='10 11'       TF='14 18' FitTwoStage 0

  DeltaT="24 28" TI='8 8' TF='18 22' FitTwoStage 1
  Alt= DeltaT="24 28 32" TI='8 8 16' TF='18 22 18' FitTwoStage 1 # Just to see Delta T=32
  Alt= DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 1
  Alt= DeltaT="24 28 32" TI='10 10 16' TF='17 21 18' FitTwoStage 1 # Just to see Delta T=32

  DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 2
  Alt= DeltaT="20 24 28" TI='10 10 10' TF='13 17 21' FitTwoStage 2
  Alt= DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 2
  Alt= DeltaT="24 28 32" TI='10 10 15' TF='17 21 17' FitTwoStage 2 # Just to see Delta T=32

  DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 3 # Preferred
  Alt= DeltaT="20 24 28" TI='10 10 10' TF='13 17 21' FitTwoStage 3
  Alt= DeltaT="16 20" TI='10 10' TF='10 13' FitTwoStage 3; echo "Ignore previous error (affects raw correlator plot only)"
  Alt= DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 3
  Alt= DeltaT="24 28 32" TI='10 10 20' TF='17 21 21' FitTwoStage 3 # Just to see Delta T=32

  DeltaT="16 20" TI='7 8' TF='10 11' FitTwoStage 4
  Alt= DeltaT="16 20 24" TI='7 8 10' TF='10 11 15' FitTwoStage 4
  Alt= DeltaT="20 24" TI='8 10' TF='11 15' FitTwoStage 4
  Alt= DeltaT="20 24 28" TI='8 10 12' TF='11 15 18' FitTwoStage 4
  Alt= DeltaT="20 24 28 32" TI='8 10 12 20' TF='11 15 18 21' FitTwoStage 4
  #DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 4
  #Alt= DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 4
  #Alt= DeltaT="20 24 28 32" TI='10 10 14 20' TF='13 17 15 21' FitTwoStage 4 # See DT 28,32

  Gamma=gXYZ
  DeltaT="20 24 28" TI='10 10 10' TF='14 18 22' FitTwoStage 1 # Preferred
  Alt= DeltaT="20 24" TI='10 10' TF='14 18' FitTwoStage 1
  Alt= DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 1
  Alt= DeltaT="24 28 32" TI='10 10 15' TF='17 21 23' FitTwoStage 1
  Alt= DeltaT="24 28 32" TI='10 10 10' TF='17 21 25' FitTwoStage 1
  Alt= DeltaT="24 28 32" TI='10 10 16' TF='17 21 18' FitTwoStage 1 # Just to see Delta T=32

  DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 2
  Alt= DeltaT="20 24 28" TI='10 10 10' TF='13 17 21' FitTwoStage 2 # Preferred
  Alt= DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 2
  Alt= DeltaT="24 28 32" TI='10 10 15' TF='17 21 23' FitTwoStage 2
  Alt= DeltaT="24 28 32" TI='10 10 10' TF='17 21 25' FitTwoStage 2 # Just to see Delta T=32

  DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 3 # Preferred
  Alt= DeltaT="20 24 28" TI='10 10 13' TF='13 17 21' FitTwoStage 3
  Alt= DeltaT="24 28" TI='10 13' TF='17 21' FitTwoStage 3
  Alt= DeltaT="24 28 32" TI='10 13 10' TF='17 21 25' FitTwoStage 3 # Just to see Delta T=32

  DeltaT="20 24" TI='10 12' TF='13 17' FitTwoStage 4 # Preferred
  Alt= DeltaT="20 24 28" TI='10 12 17' TF='13 17 21' FitTwoStage 4 # Preferred
  Alt= DeltaT="24 28" TI='12 17' TF='17 21' FitTwoStage 4
  Alt= DeltaT="24 28 32" TI='12 17 21' TF='17 21 25' FitTwoStage 4 # Just to see Delta T=32
}

############################################################

# Give me an idea of what the ratios look like

############################################################

function RatioFitsIdea()
{
  echo "M2 performing Idea ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  NumExp=2 DeltaT="16 20 24 28 32" TI='9 11 11 11 11' TF='10 13 18 22 26' FitTwoStage 0
  for Gamma in gT gXYZ; do
  NumExp=2 DeltaT="16 20 24 28 32" TI='8 10  9  8  9' TF='10 13 18 22 25' FitTwoStage 1
  NumExp=2 DeltaT="16 20 24 28 32" TI='7  8  7  8  8' TF='10 13 18 22 26' FitTwoStage 2
  NumExp=2 DeltaT="16 20 24 28 32" TI='7  7  7  8 14' TF='10 14 18 22 25' FitTwoStage 3
  NumExp=2 DeltaT="16 20 24 28 32" TI='7  7  8  7 17' TF='10 13 18 22 26' FitTwoStage 4
  done
}

############################################################

# Ratio fit ranges being tested

############################################################

function RatioFitsTest()
{
  echo "M2 performing Test ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  Gamma=gT
  UnCorr= Alt= NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 4
}

############################################################

# Simultaneously fit two-pt functions and ratios

############################################################

function RatioFitsSimul()
{
  echo "M2 performing Simul ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  Gamma=gT
  Thin=('' '' t3)
  NumExp=3
  DeltaT=(24 28 32); TI=(10 10 10); TF=(18 23 27); DoSimulFit 0
  SourcePriorFit=R3_l_h447_gT_p2_0.dt_24_28_32.Simul.corr_7_18_11_18_10_28_15_28_10_18_10_23_10_27.g5P_g5W.model.$MLUSeed.h5; SinkPriorDisp=
  NumExp=2
  #DeltaT=(24 28 32); TI=(9 9 9); TF=(17 22 26); DoSimulFit 1
  #DeltaT=(28 32); TI=(9 9); TF=(22 26); DoSimulFit 1
}

############################################################

# Input

############################################################

FitWhat="${FitWhat-R3}"
Ratio=ratioE1ZV1

qSrc=h$Heavy
qSnk=l
qSpec=s

ChooseHeavyFits

for FileSeries in ${series-disp renorm AltZV}
do
(
  case $FileSeries in
    disp)
      ChooseTwoPtFits $FileSeries
      RatioFitsBase
      [ -v UnCorr ] || UnCorr= RatioFitsBase;;

    renorm)
      ChooseTwoPtFits disp
      Ratio=ratio Renorm= NotRaw= RatioFitsBase;;

    AltZV)
      ChooseTwoPtFits disp
      Ratio=ratioAltZV Renorm= NotRaw= RatioFitsBase;;

    simul)
      RatioFitsSimul;;

    idea)
      ChooseTwoPtFits disp
      RatioFitsIdea;;

    "test")
      ChooseTwoPtFits disp
      RatioFitsTest;;

    *) echo "FileSeries $FileSeries unrecognised"; exit 1;;
  esac
)
done
