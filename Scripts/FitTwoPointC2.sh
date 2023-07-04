#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=${Ensemble:-C2}
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
  Meson=h6413_s SubDir=e2_1a TwoPointScan t=3:24:10:1,e=2 t=13:24,e=1
  p=0 Meson=s_l SubDir=e2_1a FitOptions='--Hotelling .001' TwoPointScan t=3:20:10:1,e=2 t=7:20,e=1
  p=1 Meson=s_l SubDir=e2_1a TwoPointScan t=3:22:10:1,e=2 t=7:22,e=1
  p=2 Meson=s_l SubDir=e2_1a TwoPointScan t=3:20:10:1,e=2 t=7:20,e=1
  p=3 Meson=s_l SubDir=e2_1a TwoPointScan t=3:20:10:1,e=2 t=7:20,e=1
  p=4 Meson=s_l SubDir=e2_1a TwoPointScan t=3:18:10:1,e=2 t=5:18,e=1
fi

if [ -v DoScanW ]; then
  Meson=h6413_s SubDir=e2_1b TwoPointScan t=8:24,e=2 t=10:24:10:1,e=1
  p=0 Meson=s_l SubDir=e2_1b TwoPointScan t=5:20,e=2 t=3:20:13:1,e=1
  p=1 Meson=s_l SubDir=e2_1b TwoPointScan t=6:22,e=2 t=3:22:13:1,e=1
  p=2 Meson=s_l SubDir=e2_1b TwoPointScan t=5:20,e=2 t=3:20:13:1,e=1
  p=3 Meson=s_l SubDir=e2_1b TwoPointScan t=6:20,e=2 t=3:20:13:1,e=1
  p=4 Meson=s_l SubDir=e2_1b FitOptions='--Hotelling .001' TwoPointScan t=5:18,e=2 t=3:18:13:1,e=1
fi

############################################################

declare -A ayRange
ayRange[h6413_s,0]='1.095:1.135'
ayRange[s_l,0]='0.315:0.365'
ayRange[s_l,1]='0.4:0.45'
ayRange[s_l,2]='0.45:0.56'
ayRange[s_l,3]='0.51:0.64'
ayRange[s_l,4]='0.57:0.75'
#ayRange[s_l,4]='0.48:0.8'

if [ "${DoAll+x}" = x ]; then
  # D_s
  (
  export Meson=h6413_s
  export ti='8 8'
  p=0 TI=8 TF=24 NumExp=2 TI2=13 NumExp2=1 FitTwoPoint # Preferred 27 Mar 23
  p=0 TI=8 TF=24 NumExp=2 FitTwoPoint
  FitOptions='--nopolap g5P'
  ExtraName=dispind
  p=0 TI=8 TF=24 NumExp=2 TI2=13 NumExp2=1 FitTwoPoint # Preferred 27 Mar 23
  p=0 TI=8 TF=24 NumExp=2 FitTwoPoint
  )
  # Kaon
  (
  export Meson=s_l
  export LabelTF=30
  export ti='3 3'
  p=0 TI=5 TF=20 NumExp=2 FitTwoPoint
  p=0 TI=5 TF=20 NumExp=2 TI2=6 FitOptions='--Hotelling 0' FitTwoPoint
  p=0 TI=5 TF=20 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint # Preferred
  p=0 TI=5 TF=18 NumExp=2 TI2=6 NumExp2=1 TF2=18 FitOptions='--Hotelling 0' FitTwoPoint
  p=0 TI=6 TF=18 NumExp=2 TI2=6 NumExp2=1 TF2=18 FitOptions='--Hotelling 0' FitTwoPoint
  p=1 TI=6 TF=22 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint # Preferred
  p=1 TI=6 TF=23 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint
  p=2 TI=5 TF=20 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint # Preferred
  p=2 TI=6 TF=20 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint
  p=3 TI=5 TF=20 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint
  p=3 TI=6 TF=20 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint # Preferred
  p=4 TI=5 TF=18 NumExp=2 NumExp2=1 FitTwoPoint # Preferred
  FitOptions='--guess=0.61,0.9,28,11,5500' p=4 TI=6 TF=16 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint
  )
fi

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
  aTimes=(5:20 6:22 5:20 6:20 5:18)
  aTimesW=(7:20 7:22 7:20 7:20 5:18)
  case $i in
    1) aThin=('' 2 2 2 2); aTimes[2]='5:19'; aTimes[4]='5:17';;
    2) aThin=('' 1:3:2 1:3:2 1:3:2 1:3:2); aTimes[1]='6:21'; aTimes[3]='6:19';; # Preferred
    *) unset aThin;;
  esac
  PriorFitTimes="${aTimes[0]}_${aTimesW[0]}"; PriorFitTimes="${PriorFitTimes//:/_}"
  #echo "PriorFitTimes='$PriorFitTimes'"
  yrange=0.3:0.7 SimulP s_l # Simultaneous fits of point-point data at all momenta
  (
    MultiFit="$MultiFit --nopolap g5P"
    yrange=0.3:0.7 SimulP s_l dispind
  )
  done
fi

# testing 20 Apr 2023 (covariance source)
#  export Meson=s_l
#  export LabelTF=30
#  export ti='3 3'
#p=0 TI=5 TF=20 NumExp=2 TI2=7 NumExp2=1 FitTwoPoint

set +e # Keep last (so we can source this script without failures causing exit)
