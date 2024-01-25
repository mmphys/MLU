#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=${Ensemble:-C1}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit
fi
. PlotCommon.sh
. FitTwoPoint.sh

############################################################

if [ -v DoScan ]; then
  DoScanP=
  DoScanW=
  #Meson=s_l SubDir=e2_2 TwoPointScan t=4:22:5:2,e=2 t=3:22:3:2,e=2
  #Meson=s_l SubDir=e2_1 SummaryOptions=--tablen=10 TwoPointScan t=4:22:3:2,e=2 t=6:22:3:2,e=1
fi

if [ -v DoScanP ]; then
  Meson=h6413_s SubDir=e2_1a TwoPointScan t=3:27:10:1,e=2 t=14:27,e=1
  p=0 Meson=s_l SubDir=e2_1a TwoPointScan t=3:23:10:1,e=2 t=7:23,e=1
  p=1 Meson=s_l SubDir=e2_1a TwoPointScan t=3:23:10:1,e=2 t=7:23,e=1
  p=2 Meson=s_l SubDir=e2_1a TwoPointScan t=3:20:10:1,e=2 t=5:20,e=1
  p=3 Meson=s_l SubDir=e2_1a TwoPointScan t=3:20:10:1,e=2 t=5:20,e=1
  p=4 Meson=s_l SubDir=e2_1a TwoPointScan t=3:18:10:1,e=2 t=5:18,e=1
fi

if [ -v DoScanW ]; then
  Meson=h6413_s SubDir=e2_1b TwoPointScan t=6:27,e=2 t=13:27:10:1,e=1
  p=0 Meson=s_l SubDir=e2_1b TwoPointScan t=6:23,e=2 t=3:23:13:1,e=1
  p=1 Meson=s_l SubDir=e2_1b TwoPointScan t=6:23,e=2 t=3:23:13:1,e=1
  p=2 Meson=s_l SubDir=e2_1b TwoPointScan t=5:20,e=2 t=3:20:13:1,e=1
  p=3 Meson=s_l SubDir=e2_1b TwoPointScan t=5:20,e=2 t=3:20:13:1,e=1
  p=4 Meson=s_l SubDir=e2_1b TwoPointScan t=5:18,e=2 t=3:18:13:1,e=1
fi

############################################################

declare -A ayRange
ayRange[h6413_s,0]='1.09:1.12'
ayRange[s_l,0]='0.28:0.33'
ayRange[s_l,1]='0.37:0.44'
ayRange[s_l,2]='0.42:0.52'
ayRange[s_l,3]='0.48:0.6'
ayRange[s_l,4]='0.55:0.7'
#ayRange[s_l,4]='0.48:0.8'

if [ -v DoOld ]; then
  # D_s
  (
  export Meson=h6413_s
  export ti='8 8'
  p=0 TI=6 TF=27 NumExp=2 TI2=17 TF2=27 NumExp2=1 FitTwoPoint
  )
  # Kaon
  (
  export Meson=s_l
  export LabelTF=30
  export ti='3 2'
  p=0 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 FitTwoPoint
  p=1 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 FitTwoPoint
  p=2 TI=5 TF=20 NumExp=2 TI2=5 TF2=20 NumExp2=1 FitTwoPoint
  p=3 TI=5 TF=20 NumExp=2 TI2=5 TF2=20 NumExp2=1 FitTwoPoint
  p=4 TI=5 TF=18 NumExp=2 TI2=5 TF2=18 NumExp2=1 FitTwoPoint
  )
fi

if [ -v DoAll ]; then
  # D_s
  (
  export Meson=h6413_s
  export ti='8 8'
  p=0 TI=6 TF=27 NumExp=2 TI2=14 TF2=27 NumExp2=1 FitTwoPoint # Preferred 1 Mar 23
  )
fi

# Simultaneous fits follow
(
  InDir=$PlotData/corr/2ptp2
  OutDir=$Ensemble/MELFit/2ptp2/s_l

  aFitFiles=()
  for((i=0; i < 5; ++i)); do
    aFitFiles+=($InDir/s_l_p2_${i}_g5P_g5P.fold.$DataSeed)
  done
  MultiFit="MultiFit --Hotelling 0 --overwrite --debug-signals --strict 3"
  [ -v FitOptions ] && MultiFit="$MultiFit $FitOptions"

  aTimes=(6:23 6:23 5:20 5:20 5:18)
  aTimesW=(7:23 7:23 5:20 5:20 5:18)
  PriorFitTimes="${aTimes[0]}_${aTimesW[0]}"; PriorFitTimes="${PriorFitTimes//:/_}"

if [ -v DoDisp ]; then
  yrange=0.29:0.7 SimulP s_l # Simultaneous fits of point-point data at all momenta
  yrange=0.29:0.7 ExtraName=continuum FitOptions='--dispersion continuum' SimulP s_l
fi
if [ -v DoOld ]; then
  FitEachMomPW priorPW
  FitEachMomP priorP
  aTimes=(6:23 6:23 6:20 6:20 6:18)
  aTimesW=(7:23 7:23 5:20 6:20 6:18)
  FitEachMomPW betterPW
  FitEachMomP betterP
fi
)
