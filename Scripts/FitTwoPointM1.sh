#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=${Ensemble:-M1}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit
fi
. PlotCommon.sh
. FitTwoPoint.sh

#export FitOptions="--nophat${FitOptions:+ $FitOptions}" # Enable to use p instead of p_hat

############################################################

if [ -v DoScan ]; then
  DoScanP=
  DoScanW=
fi

if [ -v DoScanP ]; then
  Meson=h${Heavy}_s SubDir=e2_1a TwoPointScan t=5:28:11:1,e=2 t=15:28,e=1
  Meson=h${Heavy}_s SubDir=P2W2a TwoPointScan t=5:28:11:1,e=2 t=9:28,e=2
  p=0 Meson=s_l SubDir=e2_1a TwoPointScan t=3:18:8:1,e=2 t=11:18,e=1
  p=0 Meson=s_l SubDir=P2W2a TwoPointScan t=3:18:8:1,e=2 t=6:18,e=2
  p=1 Meson=s_l SubDir=e2_1a TwoPointScan t=3:19:8:1,e=2 t=8:19,e=1
  p=2 Meson=s_l SubDir=e2_1a TwoPointScan t=3:19:8:1,e=2 t=9:19,e=1
  p=3 Meson=s_l SubDir=e2_1a TwoPointScan t=3:16:8:1,e=2 t=8:21,e=1
  p=4 Meson=s_l SubDir=e2_1a TwoPointScan t=4:19:7:1,e=2 t=10:21,e=1
fi

if [ -v DoScanW ]; then
  Meson=h${Heavy}_s SubDir=e2_1b TwoPointScan t=10:28,e=2 t=11:28:9:1,e=1
  Meson=h${Heavy}_s SubDir=P2W2b TwoPointScan t=10:28,e=2 t=5:28:9:1,e=2
  p=0 Meson=s_l SubDir=e2_1b TwoPointScan t=6:18,e=2 t=6:18:11:1,e=1
  p=0 Meson=s_l SubDir=P2W2b TwoPointScan t=6:18,e=2 t=3:18:8:1,e=2
  p=1 Meson=s_l SubDir=e2_1b TwoPointScan t=6:19,e=2 t=5:18:8:1,e=1
  p=2 Meson=s_l SubDir=e2_1b TwoPointScan t=6:19,e=2 t=5:18:8:1,e=1
  p=3 Meson=s_l SubDir=e2_1b TwoPointScan t=6:16,e=2 t=5:21:8:1,e=1
  p=4 Meson=s_l SubDir=e2_1b TwoPointScan t=7:19,e=2 t=6:21:10:1,e=1
fi

############################################################

declare -A ayRange
ayRange[h${Heavy}_s,0]='0.81:0.845'
ayRange[s_l,0]='0.21:0.29'
ayRange[s_l,1]='0.28:0.36'
ayRange[s_l,2]='0.32:0.43'
ayRange[s_l,3]='0.35:0.55'
ayRange[s_l,4]='0.37:0.62'

if [ "${DoAll+x}" = x ]; then
  # D_s
  (
  export Meson=h${Heavy}_s
  export ti='8 8'
  export LabelTF=30
  p=0 NumExp=2 TI=10 TF=28 TI2=15 NumExp2=1 FitTwoPoint # Best 1-exp wall fit. Preferred
  p=0 NumExp=2 TI=10 TF=28 TI2=9 FitTwoPoint # Best 2-exp wall fit
  p=0 NumExp=2 TI=10 TF=28 TI2=10 FitTwoPoint
  p=0 NumExp=2 TI=12 TF=28 TI2=10 FitTwoPoint
  # These next two 3-exp p-p fits work, but not compatible with each other
  #p=0 NumExp=3 TI=11 TF=28 TI2=10 NumExp2=2 FitTwoPoint
  #p=0 NumExp=3 TI=6 TF=28 TI2=10 NumExp2=2 FitTwoPoint
  FitOptions='--nopolap g5P'
  ExtraName=dispind
  p=0 NumExp=2 TI=10 TF=28 TI2=15 NumExp2=1 FitTwoPoint
  p=0 NumExp=2 TI=10 TF=28 TI2=10 FitTwoPoint
  )
  # Kaon
  (
  export Meson=s_l
  export LabelTF=30
  export ti='3 3'
  p=0 NumExp=2 TI=6 TF=18 TI2=11 NumExp2=1 FitTwoPoint # Preferred
  p=0 NumExp=2 TI=7 TF=18 TI2=11 NumExp2=1 FitTwoPoint # After discussion with Tobi
  p=0 NumExp=2 TI=6 TF=18 TI2=13 NumExp2=1 FitTwoPoint
  p=0 NumExp=2 TI=6 TF=18 TI2=8 NumExp2=1 FitTwoPoint
  p=0 NumExp=2 TI=6 TF=18 FitTwoPoint
  p=1 NumExp=2 TI=6 TF=19 TI2=8 NumExp2=1 FitTwoPoint
  p=1 NumExp=2 TI=7 TF=19 TI2=8 NumExp2=1 FitTwoPoint # After discussion with Tobi
  p=2 NumExp=2 TI=6 TF=19 TI2=9 NumExp2=1 FitTwoPoint
  p=2 NumExp=2 TI=7 TF=19 TI2=9 NumExp2=1 FitTwoPoint # After discussion with Tobi
  p=3 NumExp=2 TI=6 TF=16 TI2=8 TF2=21 NumExp2=1 FitTwoPoint
  p=4 NumExp=2 TI=7 TF=19 TI2=10 TF2=21 NumExp2=1 FitTwoPoint
  )
fi

# Simultaneous Kaon fits, multiple momenta
if [ -v DoNext ]; then
  InDir=$PlotData/corr/2ptp2
  OutDir=$Ensemble/MELFit/2ptp2/s_l

  aFitFiles=()
  for((i=0; i < 5; ++i)); do
    aFitFiles+=($InDir/s_l_p2_${i}_g5P_g5P.fold.$DataSeed)
  done
  MultiFit="MultiFit --Hotelling 0 --overwrite --debug-signals --strict 3"
  [ -v FitOptions ] && MultiFit="$MultiFit $FitOptions"

  # Compare unthinned with a couple of choices of thinning
  for (( i=0; i<3; ++i )); do
  aTimes=(6:18 6:19 6:19 6:16 7:19)
  aTimesW=(11:18 8:19 9:19 8:21 10:21)
  case $i in
    1) aThin=('' 2 2 2 2); aTimes[1]='6:18'; aTimes[2]='6:18';;
    2) aThin=('' 1:3:2 1:3:2 1:3:2 1:3:2); aTimes[3]='6:15'; aTimes[4]='7:18';; # Preferred
    *) unset aThin;;
  esac
  PriorFitTimes="${aTimes[0]}_${aTimesW[0]}"; PriorFitTimes="${PriorFitTimes//:/_}"
  #echo "PriorFitTimes='$PriorFitTimes'"
  (
    export yrange=0.21:0.62
    SimulP s_l # Simultaneous fits of point-point data at all momenta
    MultiFit="$MultiFit --nopolap g5P"
    SimulP s_l dispind
  )
  done
fi

set +e # Keep last (so we can source this script without failures causing exit)
