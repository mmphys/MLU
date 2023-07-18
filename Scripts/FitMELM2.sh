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

# Input

############################################################

FitWhat="${FitWhat-R3}"
Ratio=ratioE1ZV1

qSrc=h$Heavy
qSnk=l
qSpec=s

############################################################

# Two point fit choices

############################################################

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

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

aDsTIP=10
aDsTFP=26
aDsTIW=19
aDsTFW=26

UseDsWall=1
if (( UseDsWall )); then
  # Ds point-point and point-wall (though they have some tension)
  aMesonFit[h${Heavy}_s,0]=corr_${aDsTIP}_${aDsTFP}_${aDsTIW}_${aDsTFW}
  aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
else
  # Ds point-point only. Is this any better?
  aMesonFit[h${Heavy}_s,0]=corr_${aDsTIP}_${aDsTFP}
  aMesonFileOp[h${Heavy}_s,0]=g5P
fi
aMesonFileMom[h${Heavy}_s,0]=_p2_0

aKaonTIP=( 6  6  7  7  6)
aKaonTFP=(20 20 20 19 16)
aKaonTIW=( 8  7  6  6 13)
aKaonTFW=(20 20 20 19 16)

############################################################

# Choose two step fits - i.e. fit R3 using 2pt fits as input

############################################################

function ChooseTwoPtFits()
{
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
      aMesonFit[s_l,$i]=corr_6_20_6_19_7_20_7_18_6_15 # thinned 1:3:2 # Preferred
      #aMesonFit[s_l,$i]=corr_6_20_6_20_7_19_7_19_6_16 # thinned 2
      #aMesonFit[s_l,$i]=corr_6_20_6_20_7_20_7_19_6_16 # unthinned
      aMesonFileOp[s_l,$i]=g5P
      aMesonFileMom[s_l,$i]=
    done;;

    *) echo "FileSeries $FileSeries unrecognised"; exit 1;;
  esac
}

############################################################

# Two step fits - i.e. fit R3 using 2pt fits as input

############################################################

if [ -v Do2Step ]; then

# First choose which 2pt fits to use

unset FitOptions
for FileSeries in ${series-disp}
do
  echo "Performing $FitWhat fits to $Ratio for $FileSeries"
  ChooseTwoPtFits
  (
  #FitOptions='--covsrc bootstrap'
  for CorrUncorr in 0 1
  do
    # For now, these are borrowed directly from M1
    Gamma=gT
    NumExp=2 DeltaT="24 28 32" TI='11 11 11' TF='18 22 26' FitTwoStage 0
    Alt= NumExp=2 DeltaT="28 32" TI='11 11' TF='22 26' FitTwoStage 0
    Alt= NumExp=2 DeltaT="20 24 28 32" TI='10 11 11 11' TF='14 18 22 26' FitTwoStage 0
    Alt= NumExp=2 DeltaT="20 24 28   " TI='10 11 11'    TF='14 18 22' FitTwoStage 0
    Alt= NumExp=2 DeltaT="20 24      " TI='10 11'       TF='14 18' FitTwoStage 0

    NumExp=2 DeltaT="24 28" TI='8 8' TF='18 22' FitTwoStage 1
    Alt= NumExp=2 DeltaT="24 28 32" TI='8 8 16' TF='18 22 18' FitTwoStage 1 # Just to see Delta T=32
    Alt= NumExp=2 DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 1
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 16' TF='17 21 18' FitTwoStage 1 # Just to see Delta T=32

    NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 2
    Alt= NumExp=2 DeltaT="20 24 28" TI='10 10 10' TF='13 17 21' FitTwoStage 2
    Alt= NumExp=2 DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 2
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 15' TF='17 21 17' FitTwoStage 2 # Just to see Delta T=32

    NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 3 # Preferred
    Alt= NumExp=2 DeltaT="20 24 28" TI='10 10 10' TF='13 17 21' FitTwoStage 3
    Alt= NumExp=2 DeltaT="16 20" TI='10 10' TF='10 13' FitTwoStage 3
    Alt= NumExp=2 DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 3
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 20' TF='17 21 21' FitTwoStage 3 # Just to see Delta T=32

    NumExp=2 DeltaT="16 20" TI='7 8' TF='10 11' FitTwoStage 4
    Alt= NumExp=2 DeltaT="16 20 24" TI='7 8 10' TF='10 11 15' FitTwoStage 4
    Alt= NumExp=2 DeltaT="20 24" TI='8 10' TF='11 15' FitTwoStage 4
    Alt= NumExp=2 DeltaT="20 24 28" TI='8 10 12' TF='11 15 18' FitTwoStage 4
    Alt= NumExp=2 DeltaT="20 24 28 32" TI='8 10 12 20' TF='11 15 18 21' FitTwoStage 4
    #NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 4
    #Alt= NumExp=2 DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 4
    #Alt= NumExp=2 DeltaT="20 24 28 32" TI='10 10 14 20' TF='13 17 15 21' FitTwoStage 4 # See DT 28,32

    Gamma=gXYZ
    NumExp=2 DeltaT="20 24 28" TI='10 10 10' TF='14 18 22' FitTwoStage 1 # Preferred
    Alt= NumExp=2 DeltaT="20 24" TI='10 10' TF='14 18' FitTwoStage 1
    Alt= NumExp=2 DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 1
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 15' TF='17 21 23' FitTwoStage 1
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 10' TF='17 21 25' FitTwoStage 1
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 16' TF='17 21 18' FitTwoStage 1 # Just to see Delta T=32

    NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 2
    Alt= NumExp=2 DeltaT="20 24 28" TI='10 10 10' TF='13 17 21' FitTwoStage 2 # Preferred
    Alt= NumExp=2 DeltaT="24 28" TI='10 10' TF='17 21' FitTwoStage 2
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 15' TF='17 21 23' FitTwoStage 2
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 10 10' TF='17 21 25' FitTwoStage 2 # Just to see Delta T=32

    NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 3 # Preferred
    Alt= NumExp=2 DeltaT="20 24 28" TI='10 10 13' TF='13 17 21' FitTwoStage 3
    Alt= NumExp=2 DeltaT="24 28" TI='10 13' TF='17 21' FitTwoStage 3
    Alt= NumExp=2 DeltaT="24 28 32" TI='10 13 10' TF='17 21 25' FitTwoStage 3 # Just to see Delta T=32

    NumExp=2 DeltaT="20 24" TI='10 12' TF='13 17' FitTwoStage 4 # Preferred
    Alt= NumExp=2 DeltaT="20 24 28" TI='10 12 17' TF='13 17 21' FitTwoStage 4 # Preferred
    Alt= NumExp=2 DeltaT="24 28" TI='12 17' TF='17 21' FitTwoStage 4
    Alt= NumExp=2 DeltaT="24 28 32" TI='12 17 21' TF='17 21 25' FitTwoStage 4 # Just to see Delta T=32

    UnCorr=
  done
  )
done
fi

if [ -v DoIdea ]; then # Give me an idea of what the ratios look like
  FileSeries=disp; ChooseTwoPtFits; FileSeries=idea; Gamma=gXYZ
  NumExp=2 DeltaT="16 20 24 28 32" TI='9 11 11 11 11' TF='10 13 18 22 26' FitTwoStage 0
  NumExp=2 DeltaT="16 20 24 28 32" TI='8 10  9  8  9' TF='10 13 18 22 25' FitTwoStage 1
  NumExp=2 DeltaT="16 20 24 28 32" TI='7  8  7  8  8' TF='10 13 18 22 26' FitTwoStage 2
  NumExp=2 DeltaT="16 20 24 28 32" TI='7  7  7  8 14' TF='10 14 18 22 25' FitTwoStage 3
  NumExp=2 DeltaT="16 20 24 28 32" TI='7  7  8  7 17' TF='10 13 18 22 26' FitTwoStage 4
fi

if [ -v DoTest ]; then # Do manual tests here
  (
  FileSeries=disp; ChooseTwoPtFits; Gamma=gT
    UnCorr= Alt= NumExp=2 DeltaT="20 24" TI='10 10' TF='13 17' FitTwoStage 4
  )
fi

if [ -v DoSimul ]; then
  (
  echo "Performing Simultaneous $FitWhat fits to $Ratio"
  Gamma=gT
  Thin=('' '' t3)
  NumExp=3
  DeltaT=(24 28 32); TI=(10 10 10); TF=(18 23 27); DoSimulFit 0
  SourcePriorFit=R3_l_h447_gT_p2_0.dt_24_28_32.Simul.corr_7_18_11_18_10_28_15_28_10_18_10_23_10_27.g5P_g5W.model.$MLUSeed.h5; SinkPriorDisp=
  NumExp=2
  #DeltaT=(24 28 32); TI=(9 9 9); TF=(17 22 26); DoSimulFit 1
  #DeltaT=(28 32); TI=(9 9); TF=(22 26); DoSimulFit 1
  )
fi
