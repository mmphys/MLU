#!/usr/bin/env bash

export Ensemble=${Ensemble:-C1}
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
declare -A ayrangeMEL
ayrangeR3[gT,0]=0.00150:0.00190
ayrangeR3[gT,1]=0.00140:0.00165
ayrangeR3[gT,2]=0.00130:0.00160
ayrangeR3[gT,3]=0.00100:0.00160
ayrangeR3[gT,4]=0.00100:0.00150
ayrangeMEL[gT,0]=0.76:0.80
ayrangeMEL[gT,1]=0.68:0.71
ayrangeMEL[gT,2]=0.58:0.66
ayrangeMEL[gT,3]=0.50:0.64
ayrangeMEL[gT,4]=0.44:0.62
ayrangeR3[gXYZ,1]=0.00050:0.00065
ayrangeR3[gXYZ,2]=0.00038:0.00055
ayrangeR3[gXYZ,3]=0.00032:0.00046
ayrangeR3[gXYZ,4]=0.00026:0.00046
ayrangeMEL[gXYZ,1]=0.66:0.74
ayrangeMEL[gXYZ,2]=0.56:0.70
ayrangeMEL[gXYZ,3]=0.48:0.65
ayrangeMEL[gXYZ,4]=0.44:0.65

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

aMesonFit[h${Heavy}_s,0]=corr_6_27_14_27
aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
aMesonFileMom[h${Heavy}_s,0]=_p2_0

############################################################

# Choose two step fits - i.e. fit R3 using 2pt fits as input

############################################################

function ChooseTwoPtFits()
{
  case "$1" in
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
    aMesonFileMom[s_l,$i]=
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
    aMesonFileMom[s_l,$i]=_p2_$i.${Fit2ptSeries}_6_23_7_23
  done;;

  betterPW) # Simultaneous fit to PP & PW at all momenta using dispersion relation
  aMesonFit[s_l,0]=corr_6_23_7_23
  aMesonFit[s_l,1]=corr_6_23_7_23
  aMesonFit[s_l,2]=corr_6_20_5_20
  aMesonFit[s_l,3]=corr_6_20_6_20
  aMesonFit[s_l,4]=corr_6_18_6_18
  aMesonFileOp[s_l,0]=g5P_g5W
  aMesonFileMom[s_l,0]=_p2_0
  for((i = 1; i < 5; ++i)); do
    aMesonFileOp[s_l,$i]=g5P_g5W
    aMesonFileMom[s_l,$i]=_p2_$i.${Fit2ptSeries}_6_23_7_23
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
    aMesonFileMom[s_l,$i]=_p2_$i.${Fit2ptSeries}_6_23_7_23
  done;;

  betterP) # Simultaneous fit to PP at all momenta using dispersion relation
  aMesonFit[s_l,0]=corr_6_23_7_23
  aMesonFit[s_l,1]=corr_6_23
  aMesonFit[s_l,2]=corr_6_20
  aMesonFit[s_l,3]=corr_6_20
  aMesonFit[s_l,4]=corr_6_18
  aMesonFileOp[s_l,0]=g5P_g5W
  aMesonFileMom[s_l,0]=_p2_0
  for((i = 1; i < 5; ++i)); do
    aMesonFileOp[s_l,$i]=g5P_g5W
    aMesonFileMom[s_l,$i]=_p2_$i.${Fit2ptSeries}_6_23_7_23
  done;;

  *) echo "Two-point fits $1 unrecognised"; exit 1;;
esac
}

############################################################

# These were my original ratio fit ranges

############################################################

function RatioFitsBase()
{
  NumExp=2
  Gamma=gT
  DeltaT="20 24 28" TI='9 11 13' TF='14 18 18' FitTwoStage 0 # Preferred
  DeltaT="20 24 28" TI='11 14 13' TF='12 15 17' FitTwoStage 1
  DeltaT="20 24 28" TI='8 12 16' TF='12 16 19' FitTwoStage 2
  DeltaT="20 24 28" TI='8 12 16' TF='14 16 19' FitTwoStage 3
  DeltaT="20 24" TI='8 12' TF='14 16' FitTwoStage 4

  # Gamma spatial
  Gamma=gXYZ
  DeltaT="20 24 28" TI='12 15 13' TF='14 17 17' FitTwoStage 1
  DeltaT="16 20 24" TI='8 9 8' TF='11 15 19' FitTwoStage 2 # My preferred fit
  DeltaT="16 20 24" TI='8 9 8' TF='11 15 19' FitTwoStage 3 # My preferred fit
  DeltaT="16 20 24" TI='7 9 8 13' TF='11 15 19 16' FitTwoStage 4
}

############################################################

# An alternate set of ratio fit ranges ... which I never finished/used
# Seems from the following line of code I originally intended 'alt' versions of Base fit ranges?
#   FileSeries=${Fit2ptSeries}alt

############################################################

function RatioFitsAlt()
{
  Gamma=gT
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 11 13' TF='11 14 18 18' FitTwoStage 0 # Very much the same as prior
  NumExp=1 DeltaT="20 24 28" TI='11 14 13' TF='12 15 17' FitTwoStage 1 # I prefer this one
  Alt= NumExp=2 DeltaT="16 20 24 28" TI='9 11 13 13' TF='10 13 16 17' FitTwoStage 1

  Gamma=gXYZ
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 8 13' TF='11 15 19 16' FitTwoStage 2 # Just to see DT=28
  NumExp=2 DeltaT="16 20 24 28" TI='8 9 8 13' TF='11 15 19 16' FitTwoStage 3 # Just to see DT=28
}

############################################################

# A set of standardised ratio fit ranges
# Originally these seemed most useful for uncorrelated fits

############################################################

function RatioFitsStd()
{
  DeltaT="16 20 24"
  TI='8 8 8'
  TF='10 14 18'
  NumExp=2
  for (( n=0; n<5; ++n )); do
    Gamma=gT FitTwoStage $n
    if ((n)); then Gamma=gXYZ FitTwoStage $n; fi
  done
}

############################################################

# Ratio fit ranges being tested

############################################################

function RatioFitsTest()
{
  Gamma=gT
##NumExp=2 DeltaT="20 24 28" TI='10 10 13' TF='13 18 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705 FitTwoStage 1 # I prefer this one
#NumExp=2 DeltaT="20 24 28" TI='10 14 13' TF='13 15 17' yrangeR3=0.0014:0.00165 yrangeMEL=0.69:0.705 FitTwoStage 1 # I prefer this one
# This works
#NumExp=2 DeltaT="16 20" TI='6 8 12' TF='9 14 16' yrangeR3=0.0008:0.00150 yrangeMEL=0.44:0.62 FitTwoStage 4
# This works - end
}

############################################################

# Now make three-point fits

############################################################

FitWhat="R3"
Ratio=ratioE1ZV1

qSrc=h$Heavy
qSnk=l
qSpec=s

for FileSeries in ${series-old disp priorPW betterPW priorP betterP dispstd renorm}
do
(
  echo "Performing $FitWhat fits to $Ratio for $FileSeries"
  case $FileSeries in
    old | disp | priorPW | betterPW | priorP | betterP)
      ChooseTwoPtFits $FileSeries
      RatioFitsBase;;

    dispstd)
      UnCorr=
      ChooseTwoPtFits disp
      RatioFitsStd;;

    renorm)
      ChooseTwoPtFits disp
      Ratio=ratio Renorm= NotRaw= RatioFitsBase;;

    "test")
      ChooseTwoPtFits disp
      RatioFitsTest;;

    *) echo "FileSeries $FileSeries unrecognised"; exit 1;;
  esac
)
done
