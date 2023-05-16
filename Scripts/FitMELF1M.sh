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
declare -A ayrangeMEL
ayrangeR3[gT,0]=0.00068:0.000782
ayrangeR3[gT,1]=0.00060:0.00070
ayrangeR3[gT,2]=0.00056:0.00066
ayrangeR3[gT,3]=0.00054:0.000625
ayrangeR3[gT,4]=0.00050:0.00057
ayrangeR3[gT,5]=0.00045:0.0006
ayrangeR3[gT,6]=0.00044:0.00052
ayrangeMEL[gT,0]=0.527:0.537
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

declare -A aMesonFit
declare -A aMesonFileOp
declare -A aMesonFileMom

aMesonFileOp[h${Heavy}_s,0]=g5P_g5W
aMesonFileMom[h${Heavy}_s,0]=_p2_0
for((i=0;i<=MaxPSq;++i)); do
  aMesonFileOp[s_l,$i]=g5P_g5W
  aMesonFileMom[s_l,$i]=_p2_$i
done

for Fit2ptSeries in ${series-old}; do
case $Fit2ptSeries in
  old)
  aMesonFit[h385_s,0]=corr_6_29_5_29
  aMesonFit[s_l,0]=corr_10_26_7_26 # Original with Peter
  aMesonFit[s_l,1]=corr_8_28_9_28
  aMesonFit[s_l,2]=corr_8_26_8_28
  aMesonFit[s_l,3]=corr_8_28_8_28
  aMesonFit[s_l,4]=corr_10_28_8_28
  aMesonFit[s_l,5]=corr_10_28_12_28
  aMesonFit[s_l,6]=corr_10_28_12_28;;

  *) echo "Fit2ptSeries ${Fit2ptSeries} unrecognised"; exit 1;;
esac

############################################################

# Now make three-point fits

############################################################

FitWhat="R3"
Ratio=ratioE1ZV1

qSrc=h385
qSnk=l
qSpec=s
Gamma=gT

if [ -v DoAll ]; then
  FileSeries=${Fit2ptSeries}
  echo "Performing $FitWhat fits to $Ratio for $FileSeries"

  NumExp=2 DeltaT="24 28 32" TI='13 14 14' TF='16 20 24' DoFit 0

  NumExp=2 DeltaT="24 28 32" TI='13 13 13' TF='16 19 23' DoFit 1 #Alternate
  NumExp=2 DeltaT="24 28 32" TI='14 15 15' TF='17 20 24' Gamma=gXYZ DoFit 1

  NumExp=2 DeltaT="24 28 32" TI='12 13 13' TF='15 19 23' DoFit 2
  NumExp=2 DeltaT="24 28 32" TI='13 13 14' TF='16 20 24' Gamma=gXYZ DoFit 2

  NumExp=3 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' DoFit 3
  NumExp=2 DeltaT="24 28 32" TI='12 13 14' TF='16 20 24' Gamma=gXYZ DoFit 3

  NumExp=2 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' DoFit 4
  NumExp=2 DeltaT="24 28" TI='12 13 14' TF='16 20 24' Gamma=gXYZ DoFit 4

  NumExp=3 DeltaT="24 28 32" TI='11 12 13' TF='14 18 22' DoFit 5
  NumExp=2 DeltaT="24 28" TI='12 14 14' TF='16 20 24' Gamma=gXYZ DoFit 5

  NumExp=2 DeltaT="24 28 32" TI='12 15 18' TF='15 19 22' DoFit 6
  NumExp=2 DeltaT="24 28 32" TI='10 12 16' TF='16 20 24' Gamma=gXYZ DoFit 6
fi
if [ -v DoAlt ]; then
  FileSeries=${Fit2ptSeries}alt
  echo "Performing $FitWhat fits to $Ratio for $FileSeries"

  NumExp=2 DeltaT="28 32" TI='12 12' TF='20 24' DoFit 1 # Preferred
fi
if [ -v DoStd ]; then
  (
  FileSeries=${Fit2ptSeries}std
  echo "Performing $FitWhat fits to $Ratio for $FileSeries"
  #UnCorr=
  DeltaT="24 28 32"
  TI='12 12 12'
  TF='14 18 22'
  NumExp=2
  for (( n=0; n<7; ++n )); do
    DoFit $n
    if ((n)); then Gamma=gXYZ DoFit $n; fi
  done
  )
fi
done
