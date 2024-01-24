#!/usr/bin/env bash

export Ensemble=F1M
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

declare -A ayrangeR3
declare -A ayrangeMEL # Poor name. This is for the 3pt correlators - not matrix elements
ayrangeR3[gT,0]=0.00068:0.000782
ayrangeR3[gT,1]=0.00060:0.00070
ayrangeR3[gT,2]=0.00056:0.00066
ayrangeR3[gT,3]=0.00054:0.000625
ayrangeR3[gT,4]=0.00050:0.00057
ayrangeR3[gT,5]=0.00045:0.0006
ayrangeR3[gT,6]=0.00044:0.00052
ayrangeMEL[gT,0]=0.527:0.538
ayrangeMEL[gT,1]=0.486:0.496
ayrangeMEL[gT,2]=0.452:0.462
ayrangeMEL[gT,3]=0.42:0.44
ayrangeMEL[gT,4]=0.395:0.415
ayrangeMEL[gT,5]=0.36:0.4
ayrangeMEL[gT,6]=0.32:0.37
ayrangeR3[gXYZ,1]=0.000185:0.000225
ayrangeR3[gXYZ,2]=0.000155:0.000195
ayrangeR3[gXYZ,3]=0.000130:0.000170
ayrangeR3[gXYZ,4]=0.000115:0.00014
ayrangeR3[gXYZ,5]=0.000110:0.00014
ayrangeR3[gXYZ,6]=0.00008:0.00014
ayrangeMEL[gXYZ,1]=0.45:0.46
ayrangeMEL[gXYZ,2]=0.44:0.467
ayrangeMEL[gXYZ,3]=0.41:0.44
ayrangeMEL[gXYZ,4]=0.38:0.41
ayrangeMEL[gXYZ,5]=0.34:0.4
ayrangeMEL[gXYZ,6]=0.27:0.37

# y-ranges for renormalised R3 ratios
declare -A ayrangeR3R
ayrangeR3R[gT,0]=0.78:0.92
ayrangeR3R[gT,1]=0.68:0.82
ayrangeR3R[gT,2]=0.64:0.78
ayrangeR3R[gT,3]=0.58:0.72
ayrangeR3R[gT,4]=0.54:0.68
ayrangeR3R[gT,5]=0.52:0.66
ayrangeR3R[gT,6]=0.51:0.65
ayrangeR3R[gXYZ,1]=0.20:0.27
ayrangeR3R[gXYZ,2]=0.16:0.23
ayrangeR3R[gXYZ,3]=0.13:0.20
ayrangeR3R[gXYZ,4]=0.11:0.18
ayrangeR3R[gXYZ,5]=0.10:0.17
ayrangeR3R[gXYZ,6]=0.08:0.15

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
  unset FitOptionsRatio
  case "$1" in
    old) # Original fit selections
      FitOptionsRatio='C2eSrc=3'
      aMesonFit[h${Heavy}_s,0]=corr_6_29_5_29
      aMesonFit[s_l,0]=corr_10_26_7_26
      aMesonFit[s_l,1]=corr_8_28_9_28
      aMesonFit[s_l,2]=corr_8_26_8_28
      aMesonFit[s_l,3]=corr_8_28_8_28
      aMesonFit[s_l,4]=corr_10_28_8_28
      aMesonFit[s_l,5]=corr_10_28_12_28
      aMesonFit[s_l,6]=corr_10_28_12_28
      for((i=0;i<=MaxPSq;++i)); do
        aMesonFileOp[s_l,$i]=g5P_g5W
        aMesonFileMom[s_l,$i]=_p2_$i
      done;;

    better) # July '23 fit selections
      aMesonFit[h${Heavy}_s,0]=corr_13_28_13_29
      aMesonFit[s_l,0]=corr_7_31_8_20
      aMesonFit[s_l,1]=corr_10_26_9_21
      aMesonFit[s_l,2]=corr_9_21_8_21
      aMesonFit[s_l,3]=corr_9_21_8_21
      aMesonFit[s_l,4]=corr_9_19_10_17
      aMesonFit[s_l,5]=corr_7_28_12_28
      aMesonFit[s_l,6]=corr_7_23_12_20
      for((i=0;i<=MaxPSq;++i)); do
        aMesonFileOp[s_l,$i]=g5P_g5W
        aMesonFileMom[s_l,$i]=_p2_$i
      done;;

    disp) # Simultaneous fit to PP at all momenta using dispersion relation
      aMesonFit[h${Heavy}_s,0]=corr_13_28_13_29
      for((i=0;i<=MaxPSq;++i)); do
        aMesonFit[s_l,$i]=corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23
        aMesonFileOp[s_l,$i]=g5P
        aMesonFileMom[s_l,$i]=
      done;;

    *) echo "Two-point fits $1 unrecognised"; exit 1;;
  esac
}

############################################################

# These were my original ratio fit ranges

############################################################

function RatioFitsOld()
{
  echo "F1M performing Old ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  Gamma=gT
  NumExp=2 DeltaT="24 28 32" TI='13 14 14' TF='16 20 24' FitTwoStage 0

  NumExp=2 DeltaT="24 28 32" TI='13 13 13' TF='16 19 23' FitTwoStage 1 #Alternate
  NumExp=2 DeltaT="24 28 32" TI='14 15 15' TF='17 20 24' Gamma=gXYZ FitTwoStage 1

  NumExp=2 DeltaT="24 28 32" TI='12 13 13' TF='15 19 23' FitTwoStage 2
  NumExp=2 DeltaT="24 28 32" TI='13 13 14' TF='16 20 24' Gamma=gXYZ FitTwoStage 2

  NumExp=3 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' FitTwoStage 3
  NumExp=2 DeltaT="24 28 32" TI='12 13 14' TF='16 20 24' Gamma=gXYZ FitTwoStage 3

  NumExp=2 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' FitTwoStage 4
  NumExp=2 DeltaT="24 28" TI='12 13 14' TF='16 20 24' Gamma=gXYZ FitTwoStage 4

  NumExp=3 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' FitTwoStage 5
  NumExp=2 DeltaT="24 28" TI='12 14 14' TF='16 20 24' Gamma=gXYZ FitTwoStage 5

  NumExp=2 DeltaT="24 28 32" TI='12 15 18' TF='15 19 22' FitTwoStage 6
  NumExp=2 DeltaT="24 28 32" TI='10 12 16' TF='16 20 24' Gamma=gXYZ FitTwoStage 6
}

############################################################

# These were my original ratio fit ranges for uncorrelated

############################################################

function RatioFitsStd()
{
  UnCorr=
  echo "F1M performing Std ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  DeltaT="24 28 32"
  TI='12 12 12'
  TF='14 18 22'
  NumExp=2
  for (( n=0; n<=MaxPSq; ++n )); do
    Gamma=gT FitTwoStage $n
    ((n)) && Gamma=gXYZ FitTwoStage $n
  done
}

############################################################

# These are my latest ratio fit ranges
# These were chosen assuming we're fitting to the renormalised ratios, i.e.
#   Ratio=ratio; Renorm=

############################################################

function RatioFitsDisp()
{
  echo "F1M performing Disp ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  Gamma=gT
  NumExp=3
  DeltaT="24 28" TI='13 13' TF='16 20' FitTwoStage 0
  Alt= DeltaT="20 24 28" TI='12 13 13' TF='13 16 20' FitTwoStage 0
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 0
  Alt= DeltaT="28 32" TI='14 15' TF='19 22' FitTwoStage 0
  Alt= DeltaT="24 28 32" TI='14 14 15' TF='15 19 22' FitTwoStage 0

  DeltaT="24 28 32" TI='13 13 13' TF='16 20 24' FitTwoStage 1
  Alt= DeltaT="24 28" TI='13 13' TF='16 20' FitTwoStage 1
  Alt= DeltaT="28 32" TI='13 13' TF='20 24' FitTwoStage 1

  DeltaT="24 28" TI='13 13' TF='16 20' FitTwoStage 2
  Alt= DeltaT="20 24 28" TI='12 13 13' TF='13 16 20' FitTwoStage 2
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 2

  DeltaT="20 24 28" TI='12 13 13' TF='13 16 20' FitTwoStage 3
  Alt= DeltaT="24 28" TI='13 13' TF='16 20' FitTwoStage 3
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 3
  Alt= DeltaT="28 32" TI='13 13' TF='20 24' FitTwoStage 3

  DeltaT="24 28" TI='13 13' TF='16 20' FitTwoStage 4
  Alt= DeltaT="28 32" TI='13 13' TF='20 24' FitTwoStage 4
  Alt= DeltaT="24 28 32" TI='13 13 13' TF='16 20 24' FitTwoStage 4
  Alt= DeltaT="20 24 28 32" TI='12 13 13 13' TF='13 16 20 24' FitTwoStage 4

  DeltaT="24 28" TI='13 13' TF='16 20' FitTwoStage 5
  Alt= DeltaT="24 28 32" TI='13 13 13' TF='16 20 24' FitTwoStage 5
  Alt= DeltaT="20 24 28" TI='12 13 13' TF='13 16 20' FitTwoStage 5
  Alt= DeltaT="20 24 28 32" TI='12 13 13 13' TF='13 16 20 24' FitTwoStage 5

  DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 6
  Alt= DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 6
  Alt= DeltaT="24 28 32" TI='13 16 20' TF='16 20 24' FitTwoStage 6

  Gamma=gXYZ
  DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 1
  Alt= DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 1
  Alt= DeltaT="24 28 32" TI='13 16 20' TF='16 20 24' FitTwoStage 1
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 1
  Alt= DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 1
  Alt= DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 1

  DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 2
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 2
  Alt= DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 2
  Alt= DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 2

  DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 3
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 3
  Alt= DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 3
  Alt= DeltaT="24 28 32" TI='13 16 20' TF='16 20 24' FitTwoStage 3
  Alt= DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 3
  Alt= DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 3

  DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 4
  Alt= DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 4
  Alt= DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 4
  Alt= DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 4
  Alt= DeltaT="24 28 32" TI='13 16 20' TF='16 20 24' FitTwoStage 4
  Alt= DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 4

  DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 5
  Alt= DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 5
  Alt= DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 5
  Alt= DeltaT="24 28 32" TI='13 16 20' TF='16 20 24' FitTwoStage 5
  Alt= DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 5

  DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 6
  Alt= DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 6
  Alt= DeltaT="24 28 32" TI='13 16 20' TF='16 20 24' FitTwoStage 6
  Alt= DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 6
}

############################################################

# Ratio fit ranges being tested

############################################################

function RatioFitsTest()
{
  echo "F1M performing Test ${UnCorr+un}corr $FitWhat fits to $Ratio for $FileSeries"
  Gamma=gT
  NumExp=3
  #Ratio=ratio; Renorm=; NotRaw=
  #Ratio=ratioZV1; NotRaw=; yrangeR3='0.85:0.95'
  #Ratio=ratioE1; yrangeR3='0.00055:0.00061'
  #DeltaT=20 TI=10 TF=14 FitTwoStage 1
  #DeltaT=24 TI=12 TF=16 FitTwoStage 1
  #DeltaT=28 TI=13 TF=20 FitTwoStage 1
  #DeltaT=32 TI=13 TF=24 FitTwoStage 1
  #DeltaT='20 24' TI='12 13' TF='13 16' FitTwoStage 1
  #DeltaT="20 24 28" TI='12 13 16' TF='13 16 20' FitTwoStage 1
  #DeltaT="20 24 28 32" TI='12 13 16 20' TF='13 16 20 24' FitTwoStage 1
  DeltaT="24 28 32" TI='13 13 13' TF='16 20 24' FitTwoStage 1
  #DeltaT="24 28" TI='13 16' TF='16 20' FitTwoStage 1
  #DeltaT="28 32" TI='16 20' TF='20 24' FitTwoStage 1
}

############################################################

# Now make three-point fits

############################################################

FitWhat="R3"
Ratio=ratioE1ZV1

qSrc=h385
qSnk=l
qSpec=s

for FileSeries in ${series-old better oldstd renorm AltZV Jan24}
do
(
  case $FileSeries in
    old)
      ChooseTwoPtFits old
      RatioFitsOld;;
 
    oldstd)
      ChooseTwoPtFits old
      RatioFitsStd;;
 
    # Not sure either of better/disp are any good
    better | disp)
      ChooseTwoPtFits "$FileSeries"
      RatioFitsDisp;;

    renorm)
      ChooseTwoPtFits disp
      Ratio=ratio Renorm= NotRaw= RatioFitsDisp;;

    AltZV | Jan24)
      ChooseTwoPtFits disp
      Ratio=ratio$FileSeries Renorm= NotRaw= RatioFitsDisp;;

    "test")
      ChooseTwoPtFits disp
      RatioFitsTest;;

    *) echo "FileSeries $FileSeries unrecognised"; exit 1;;
  esac
)
done
