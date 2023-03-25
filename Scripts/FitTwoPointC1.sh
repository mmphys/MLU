#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=${Ensemble:-C1}

declare -A ayRange
ayRange[h6413_s,0]='1.08:1.12'
ayRange[s_l,0]='0.28:0.33'
ayRange[s_l,1]='0.37:0.44'
ayRange[s_l,2]='0.42:0.52'
ayRange[s_l,3]='0.48:0.6'
ayRange[s_l,4]='0.55:0.7'
#ayRange[s_l,4]='0.48:0.8'

#export FitOptions="--nophat" # Enable to use p instead of p_hat

############################################################

# Execute command

# Input
#  Cmd:     Command to execute
#  LogFile: Where to redirect output

############################################################

function DoCmd()
{
  #echo $Cmd
  echo $Cmd &>> $LogFile
  if ! eval $Cmd &>> $LogFile; then
    echo "Error ${PIPESTATUS[0]} executing $Cmd"
  fi
}

############################################################

# Simultaneous fits of point-point data at all momenta

############################################################

function SimulP()
{
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  sTimes="${aTimes[*]}"
  sTimes=${sTimes// /_}
  sTimes=${sTimes//:/_}
  Cmd="$MultiFit -e 2 -N 24 -o $OutDir/s_l"
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
  GetColumnValues $BaseFile.$Seed.h5 "E_0(n^2=0)=" $EKeys

  # Now plot it
  for((i=0; i < ${#aFitFiles[@]}; ++i))
  do
    ThisTF=${aTimes[i]}; ThisTF=${ThisTF##*:}; ThisTF=$((ThisTF+3))
    tf="$tf${tf:+ }$ThisTF"
    title="$title${title:+ }n^2=$i"
    files="$files${files:+ }${aFitFiles[i]}.txt"
  done
  Cmd="tf='$tf' title='$title' files='$files'"
  [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
  Cmd="$Cmd save='$BaseFile.$Seed' yrange=0.3:0.7 plottd.sh ${BaseFile}_td.$Seed.txt"
  DoCmd
  for((i=0; i < ${#aFitFiles[@]}; ++i))
  do
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
}

############################################################

# For each non-zero momentum
# fit pp & pw data using E(n^2=0) and lattice dispersion relation

# Input:
#  aTimes:  array of point-point start:stop fit times for each momenta (0'th ignored)
#  aTimesW: array of point-wall  start:stop fit times for each momenta (0'th ignored)
#  PriorFitTimes: fit to use for zero momentum

############################################################

function FitEachMomPW()
{
  local OutGroup=$1
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local ZeroMomFit="$OutDir/s_l_p2_0.corr_$PriorFitTimes.g5P_g5W.model.1835672416.h5"
  for((i=1; i < 5; ++i))
  do
    sTimes="${aTimes[i]} ${aTimesW[i]}"
    sTimes=${sTimes// /_}
    sTimes=${sTimes//:/_}
    BaseFile="$OutDir/s_l_p2_$i.${OutGroup}_$PriorFitTimes"
    Cmd="$MultiFit -e 2 -N 24 -o $BaseFile $ZeroMomFit"
    ThisFile="${aFitFiles[i]}.h5"
    Cmd="$Cmd $ThisFile,t=${aTimes[i]}"
    Cmd="$Cmd ${ThisFile//g5P_g5P/g5P_g5W},t=${aTimesW[i]},e=1"
    BaseFile="$BaseFile.corr_$sTimes.g5P_g5W.model"
    LogFile="$BaseFile.$Seed.log"
    DoCmd

    GetColumnValues $BaseFile.$Seed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    ThisFile="${aFitFiles[i]}.txt"
    Cmd="title='point-point point-wall' files='$ThisFile ${ThisFile//g5P_g5P/g5P_g5W}'"
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd save='$BaseFile.$Seed'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$Seed.txt"
    DoCmd
  done
}

############################################################

# For each non-zero momentum
# fit pp & pw data using E(n^2=0) and lattice dispersion relation

# Input:
#  aTimes:  array of point-point start:stop fit times for each momenta (0'th ignored)
#  PriorFitTimes: fit to use for zero momentum

############################################################

function FitEachMomP()
{
  local OutGroup=$1
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local ZeroMomFit="$OutDir/s_l_p2_0.corr_$PriorFitTimes.g5P_g5W.model.1835672416.h5"
  for((i=1; i < 5; ++i))
  do
    sTimes="${aTimes[i]}"
    sTimes=${sTimes//:/_}
    BaseFile="$OutDir/s_l_p2_$i.${OutGroup}_$PriorFitTimes"
    Cmd="$MultiFit -e 2 -N 24 -o $BaseFile $ZeroMomFit"
    Cmd="$Cmd ${aFitFiles[i]}.h5,t=${aTimes[i]}"
    BaseFile="$BaseFile.corr_$sTimes.g5P_g5W.model"
    LogFile="$BaseFile.$Seed.log"
    DoCmd

    GetColumnValues $BaseFile.$Seed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    Cmd="title='point-point n^2=$i' files='${aFitFiles[i]}.txt'"
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd save='$BaseFile.$Seed'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$Seed.txt"
    DoCmd
  done
}

############################################################

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
  (
  export Meson=h6413_s
  export ti='8 8'
  p=0 TI=6 TF=27 NumExp=2 TI2=14 TF2=27 NumExp2=1 yrange="${ayRange[$Meson,0]}" FitTwoPoint.sh # Preferred 1 Mar 23
  p=0 TI=6 TF=27 NumExp=2 TI2=17 TF2=27 NumExp2=1 FitTwoPoint.sh
  )
  # Kaon
  (
  export Meson=s_l
  export LabelTF=30
  export ti='3 2'
  p=0 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 yrange="${ayRange[$Meson,0]}" FitTwoPoint.sh
  p=1 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 yrange="${ayRange[$Meson,1]}" FitTwoPoint.sh
  p=2 TI=5 TF=20 NumExp=2 TI2=5 TF2=20 NumExp2=1 yrange="${ayRange[$Meson,2]}" FitTwoPoint.sh
  p=3 TI=5 TF=20 NumExp=2 TI2=5 TF2=20 NumExp2=1 yrange="${ayRange[$Meson,3]}" FitTwoPoint.sh
  p=4 TI=5 TF=18 NumExp=2 TI2=5 TF2=18 NumExp2=1 yrange="${ayRange[$Meson,4]}" FitTwoPoint.sh
  )
fi

. PlotCommon.sh
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit
fi

mkdir -p $Ensemble/MELFit
cd $Ensemble/MELFit

InDir=$PlotData/corr/2ptp2
OutDir=2ptp2/s_l

aFitFiles=()
for((i=0; i < 5; ++i)); do
  aFitFiles+=($InDir/s_l_p2_${i}_g5P_g5P.fold.$Seed)
done
MultiFit="MultiFit --Hotelling 0 --overwrite --debug-signals --strict 3"
[ -v FitOptions ] && MultiFit="$MultiFit $FitOptions"

if [ "${DoAll+x}" = x ]; then
  aTimes=(6:23 6:23 5:20 5:20 5:18)
  aTimesW=(7:23 7:23 5:20 5:20 5:18)
  PriorFitTimes="6_23_7_23"

  SimulP # Simultaneous fits of point-point data at all momenta
  FitEachMomPW priorPW
  FitEachMomP priorP
  aTimes=(6:23 6:23 6:20 6:20 6:18)
  aTimesW=(7:23 7:23 5:20 6:20 6:18)
  FitEachMomPW betterPW
  FitEachMomP betterP

fi
