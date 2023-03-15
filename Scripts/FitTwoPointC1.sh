#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=C1

if [ -v DoScan ]; then
  DoScanP=
  DoScanW=
  #Meson=s_l SubDir=e2_2 TwoPointScan.sh t=4:22:5:2,e=2 t=3:22:3:2,e=2
  #Meson=s_l SubDir=e2_1 SummaryOptions=--tablen=10 TwoPointScan.sh t=4:22:3:2,e=2 t=6:22:3:2,e=1
fi

if [ -v DoScanP ]; then
  Meson=h6413_s SubDir=e2_1a TwoPointScan.sh t=3:27:10:1,e=2 t=14:27,e=1
  p=0 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:23:10:1,e=2 t=7:23,e=1
  p=1 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:23:10:1,e=2 t=7:23,e=1
  p=2 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:20:10:1,e=2 t=5:20,e=1
  p=3 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:20:10:1,e=2 t=5:20,e=1
  p=4 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:18:10:1,e=2 t=5:18,e=1
fi

if [ -v DoScanW ]; then
  Meson=h6413_s SubDir=e2_1b TwoPointScan.sh t=6:27,e=2 t=13:27:10:1,e=1
  p=0 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=6:23,e=2 t=3:23:13:1,e=1
  p=1 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=6:23,e=2 t=3:23:13:1,e=1
  p=2 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=5:20,e=2 t=3:20:13:1,e=1
  p=3 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=5:20,e=2 t=3:20:13:1,e=1
  p=4 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=5:18,e=2 t=3:18:13:1,e=1
fi

if [ "${DoAll+x}" = x ]; then
  # D_s
  p=0 TI=6 TF=27 NumExp=2 TI2=14 TF2=27 NumExp2=1 Meson=h6413_s yrange='1.08:1.12' ti='8 8' FitTwoPoint.sh # Preferred 1 Mar 23
  p=0 TI=6 TF=27 NumExp=2 TI2=17 TF2=27 NumExp2=1 Meson=h6413_s yrange='1.08:1.12' ti='8 8' FitTwoPoint.sh

  # Kaon
  export LabelTF=30
  export ti='3 2'
  p=0 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 Meson=s_l yrange='0.28:0.33' FitTwoPoint.sh
  p=1 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 Meson=s_l yrange='0.37:0.44' FitTwoPoint.sh
  p=2 TI=5 TF=20 NumExp=2 TI2=5 TF2=20 NumExp2=1 Meson=s_l yrange='0.42:0.52' FitTwoPoint.sh
  p=3 TI=5 TF=20 NumExp=2 TI2=5 TF2=20 NumExp2=1 Meson=s_l yrange='0.48:0.6' FitTwoPoint.sh
  p=4 TI=5 TF=18 NumExp=2 TI2=5 TF2=18 NumExp2=1 Meson=s_l yrange='0.48:0.8' FitTwoPoint.sh
  unset ti
fi

function DoCmd()
{
  #echo $Cmd
  echo $Cmd &>> $LogFile
  eval $Cmd &>> $LogFile
}

. PlotCommon.sh
if [ -d $Ensemble ]; then
mkdir -p $Ensemble/MELFit
cd $Ensemble/MELFit

InDir=$PlotData/corr/2ptp2
OutDir=2ptp2/s_l

aFitFiles=()
for((i=0; i < 5; ++i)); do
  aFitFiles+=($InDir/s_l_p2_${i}_g5P_g5P.fold.$Seed)
done
aTimes=(6:23 6:23 5:20 5:20 5:18)
sTimes="${aTimes[*]}"
sTimes=${sTimes// /_}
sTimes=${sTimes//:/_}

# Perform Fit
Cmd="MultiFit -e 2 -N 24 --Hotelling 0 --overwrite --debug-signals --strict 3 -o $OutDir/s_l"
for((i=0; i < ${#aFitFiles[@]}; ++i)); do
  Cmd="$Cmd ${aFitFiles[i]}.h5,t=${aTimes[i]}"
done
BaseFile="$OutDir/s_l.corr_$sTimes.g5P.model"
LogFile="$BaseFile.$Seed.log"
DoCmd

# Get the energy difference
for((i=0; i < ${#aFitFiles[@]}; ++i)); do
  EKeys=$EKeys${EKeys:+,}s_l_p2_${i}-E0
done
if ColumnValues=$(GetColumn --exact ChiSqPerDof,pValueH,$EKeys $BaseFile.$Seed.h5)
#if ColumnValues=$(GetColumn --exact ChiSqPerDof,pValueH --partial E0 $BaseFile.$Seed.h5)
then
  #echo "OK: $ColumnValues"
  ColumnValues=($ColumnValues)
  E0="${ColumnValues[@]:16:8}"
  RefText="E_0(n^2=0)=${ColumnValues[@]:17:1}, χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})"
else
  LastError=${PIPESTATUS[0]}
  echo "Error $LastError: $ColumnValues"
  unset RefText
fi

# Now plot it
for((i=0; i < ${#aFitFiles[@]}; ++i)); do
  ThisTF=${aTimes[i]}; ThisTF=${ThisTF##*:}; ThisTF=$((ThisTF+3))
  tf="$tf${tf:+ }$ThisTF"
  title="$title${title:+ }n^2=$i"
  files="$files${files:+ }${aFitFiles[i]}.txt"
done
Cmd="tf='$tf' title='$title' files='$files'"
[ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
Cmd="$Cmd save='$BaseFile.$Seed' yrange=0.3:0.7 plottd.sh ${BaseFile}_td.$Seed.txt"
DoCmd
for((i=0; i < ${#aFitFiles[@]}; ++i)); do
  tf=${aTimes[i]}; tf=${tf##*:}; tf=$((tf+3))
  title="n^2=$i"
  files="${aFitFiles[i]}.txt"
  Cmd="title='$title' files='$files'"
  #Cmd="$Cmd tf='$tf'"
  [ -v RefText ] && Cmd="$Cmd RefText='E_0(n^2=${i})=${ColumnValues[@]:$((17+i*8)):1}, χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})'"
  Cmd="$Cmd save='$OutDir/s_l_p2_${i}.corr_$sTimes.g5P.model.$Seed'"
  Cmd="$Cmd mmin=$i mmax=$i plottd.sh ${BaseFile}_td.$Seed.txt"
  DoCmd
done

fi
