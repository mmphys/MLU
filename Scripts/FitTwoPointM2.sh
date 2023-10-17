#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=${Ensemble:-M2}
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
  Meson=h${Heavy}_s SubDir=e2_1a TwoPointScan t=8:25:6:3,e=2 t=19:26,e=1
  Meson=h${Heavy}_s SubDir=P2W2a TwoPointScan t=8:25:6:3,e=2 t=12:26,e=2
  p=0 Meson=s_l SubDir=e2_1a TwoPointScan t=4:19:6:3,e=2 t=8:20,e=1
  p=1 Meson=s_l SubDir=e2_1a TwoPointScan t=4:19:6:3,e=2 t=7:20,e=1
  p=2 Meson=s_l SubDir=e2_1a TwoPointScan t=4:19:6:3,e=2 t=6:20,e=1
  p=3 Meson=s_l SubDir=e2_1a TwoPointScan t=4:18:6:3,e=2 t=6:19,e=1
  p=4 Meson=s_l SubDir=e2_1a TwoPointScan t=4:15:6:3,e=2 t=13:16,e=1
  p=4 Meson=s_l SubDir=P2W2a TwoPointScan t=4:15:6:3,e=2 t=8:12,e=1
fi

if [ -v DoScanW ]; then
  Meson=h${Heavy}_s SubDir=e2_1b TwoPointScan t=10:26,e=2 t=18:25:3:3,e=1
  Meson=h${Heavy}_s SubDir=P2W2b TwoPointScan t=10:26,e=2 t=11:25:3:3,e=2
  p=0 Meson=s_l SubDir=e2_1b TwoPointScan t=6:20,e=2 t=7:19:3:3,e=1
  p=1 Meson=s_l SubDir=e2_1b TwoPointScan t=6:20,e=2 t=6:19:3:3,e=1
  p=2 Meson=s_l SubDir=e2_1b TwoPointScan t=7:20,e=2 t=5:19:3:3,e=1
  p=3 Meson=s_l SubDir=e2_1b TwoPointScan t=7:19,e=2 t=5:18:3:3,e=1
  p=4 Meson=s_l SubDir=e2_1b TwoPointScan t=6:16,e=2 t=12:15:3:3,e=1
  p=4 Meson=s_l SubDir=P2W2b TwoPointScan t=6:16,e=2 t=7:11:3:3,e=1
fi

############################################################

# y axis ranges for 2pt fit plots
declare -A ayRange
ayRange[h${Heavy}_s,0]='0.815:0.845'
ayRange[s_l,0]='0.22:0.28'
ayRange[s_l,1]='0.28:0.36'
ayRange[s_l,2]='0.32:0.43'
ayRange[s_l,3]='0.34:0.55'
ayRange[s_l,4]='0.37:0.62'

if [ "${DoOld+x}" = x ]; then
  # D_s
  (
  export Meson=h${Heavy}_s
  export NumExp=2
  export LabelTF=30
  export ti='8 8'
  p=0 TI=10 TF=26 TI2=12 FitTwoPoint # Best 2-exp wall fit - but still a bit sh*t
  p=0 TI=10 TF=26 TF2=0 FitTwoPoint # Point-point only
  )
  # Kaon
  (
  export Meson=s_l
  export NumExp=2
  export LabelTF=30
  export ti='3 3'
  p=0 TI=6 TF=20 TI2=8 NumExp2=1 FitTwoPoint # Preferred
  p=1 TI=6 TF=20 TI2=7 NumExp2=1 FitTwoPoint
  p=2 TI=7 TF=20 TI2=6 NumExp2=1 FitTwoPoint # Preferred
  p=2 TI=6 TF=20 TI2=6 NumExp2=1 FitTwoPoint # First attempt - dispersion relation not good
  p=3 TI=7 TF=19 TI2=6 NumExp2=1 FitTwoPoint # Preferred
  p=3 TI=6 TF=19 TI2=6 NumExp2=1 FitTwoPoint # First attempt - dispersion relation not good
  p=4 TI=6 TF=16 TI2=13 NumExp2=1 FitTwoPoint
  p=4 TI=6 TF=16 TI2=8 TF2=12 NumExp2=1 FitTwoPoint
  )
fi

if [ -v DoTest ]; then
  (
  export Meson=s_l
  export NumExp=2
  export LabelTF=30
  export ti='3 3'
  p=3 TI=7 TF=19 TI2=9 NumExp2=1 FitTwoPoint # Preferred
  )
fi

if [ -v DoAll ]; then
  # D_s
  (
  export Meson=h${Heavy}_s
  export NumExp=2
  export LabelTF=30
  export ti='8 8'
  p=0 TI=10 TF=26 TI2=19 NumExp2=1 FitTwoPoint # Best 1-exp wall fit. Preferred
  )
fi

# Simultaneous Kaon fits, multiple momenta
if [ -v DoDisp ]; then
  InDir=$PlotData/corr/2ptp2
  OutDir=$Ensemble/MELFit/2ptp2/s_l

  aFitFiles=()
  for((i=0; i < 5; ++i)); do
    aFitFiles+=($InDir/s_l_p2_${i}_g5P_g5P.fold.$DataSeed)
  done
  MultiFit="MultiFit --Hotelling 0 --overwrite --debug-signals --strict 3"
  [ -v FitOptions ] && MultiFit="$MultiFit $FitOptions"

  # Compare unthinned with a couple of choices of thinning
  aTimes=(6:20 6:20 7:20 7:19 6:16)
  for (( i=0; i<3; ++i )); do
  (
    case $i in
      # Preferred thinning: 1:3:2 on non-zero momenta
      2) aThin=('' 1:3:2 1:3:2 1:3:2 1:3:2);; # aTimes=(6:20 6:19 7:20 7:18 6:15);;
      # Others
      1) aThin=('' 2 2 2 2);; # aTimes[3]='6:18';;
      *) unset aThin;;
    esac
    yrange=0.22:0.55 SimulP s_l # Simultaneous fits of point-point data at all momenta
  )
  done
fi

set +e # Keep last (so we can source this script without failures causing exit)
