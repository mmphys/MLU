#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=${Ensemble:-F1M}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit
fi
. PlotCommon.sh
. FitTwoPoint.sh

############################################################

if [ -v DoScan ]; then
(
  p=0 Meson=h${Heavy}_s SubDir=P3W3 NumExp=3 NumExp2=3 StdTwoPointScan 6 29 5 29
  SubDir=e2_1
  p=0 Meson=h${Heavy}_s NumExp2=2 StdTwoPointScan 13 28 13 29
  Meson=s_l
  p=0 StdTwoPointScan 7 31 8 20
  p=0 SubDir=alt1 StdTwoPointScan 8 31 8 20
  p=1 StdTwoPointScan 10 26 9 21
  p=1 SubDir=alt1 StdTwoPointScan 9 26 9 21
  p=2 StdTwoPointScan 9 21 8 21
  p=2 SubDir=alt1 StdTwoPointScan 10 21 8 21
  p=3 StdTwoPointScan 9 21 8 21
  p=3 SubDir=alt1 StdTwoPointScan 9 21 12 21
  p=4 StdTwoPointScan 7 19 10 17
  p=4 SubDir=alt1 StdTwoPointScan 9 19 12 17
  p=5 StdTwoPointScan 7 28 12 28
  p=6 StdTwoPointScan 7 23 12 20
)
fi

(
  Meson=s_l
  SubDir=e2_1
)

############################################################

# y axis ranges for 2pt fit plots
declare -A ayRange
ayRange[h385_s,0]='0.72:0.739'
#ayRange[s_l,0]='0.187:0.197' # Zoom in
ayRange[s_l,0]='0.187:0.212' # Zoom out
ayRange[s_l,1]='0.225:0.26'
ayRange[s_l,2]='0.25:0.30'
ayRange[s_l,3]='0.28:0.32'
ayRange[s_l,4]='0.3:0.385'
ayRange[s_l,5]='0.32:0.40'
ayRange[s_l,6]='0.29:0.45'

LabelTF=40
LabelTI=5
NumExp=2

if [ -v DoOld ]; then
  DoAxial=
  # My original (not so good) fit choices
  # D_s
  p=0 TI=6 TF=29 NumExp=3 TI2=5 Meson=h385_s FitTwoPoint

  # Kaon
  p=0 TI=10 TF=26 TI2=7 Meson=s_l LabelTI=6 FitTwoPoint
  p=1 TI=8 TF=28 TI2=9 NumExp2=1 Meson=s_l FitTwoPoint
  p=2 TI=8 TF=26 TI2=8 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
  p=3 TI=8 TF=28 TI2=8 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
  p=3 TI=8 TF=28 TI2=6 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
  #p=3 TI=8 TF=28 TI2=6 TF2=22 NumExp2=1 Meson=s_l FitTwoPoint #also good
  p=4 TI=10 TF=28 TI2=8 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
  #p=4 TI=10 TF=25 TI2=8 TF2=23 NumExp2=1 Meson=s_l FitTwoPoint #also good
  p=5 TI=10 TF=28 TI2=12 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
  p=6 TI=10 TF=28 TI2=12 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint

  ###############################
  # Better fit choices, July 2023 - these are not as good as the one below
  # D_s
  #p=0 TI=14 TF=28 TF2=29 Meson=h385_s FitTwoPoint # Also good
  #p=0 TI=16 TF=40 TF2=29 Meson=h385_s FitTwoPoint #Shite

  # Kaon
  (
  Meson=s_l
  NumExp2=1
  p=0 TI=7 TF=31 TI2=8 TF2=20 FitTwoPoint # Preferred
  p=0 TI=8 TF=31 TF2=20 FitTwoPoint
  p=1 TI=10 TF=26 TI2=9 TF2=21 FitTwoPoint # Preferred
  p=1 TI=9  TF=26 TI2=9 TF2=21 FitTwoPoint
  p=2 TI=9  TF=21 TI2=8 FitTwoPoint # Preferred
  p=2 TI=10 TF=21 TI2=8 FitTwoPoint
  p=3 TI=9 TF=21 TI2=8 FitTwoPoint # Preferred
  p=3 TI=9 TF=21 TI2=12 FitTwoPoint
  p=4 TI=9 TF=19 TI2=10 TF2=17 FitTwoPoint # Preferred
  p=4 TI=9 TF=19 TI2=12 TF2=17 FitTwoPoint
  p=5 TI=7 TF=28 TI2=12 FitTwoPoint
  p=6 TI=7 TF=23 TI2=12 TF2=20 LabelTF=32 FitTwoPoint
  )
fi

# Testing
if [ -v DoTest ]; then
(
  Meson=s_l
  NumExp2=1
  LabelTI=6
)
fi

if [ -v DoAll ]; then
  ###############################
  # Better fit choices, July 2023
  # D_s
  LabelTI= LabelTF=35 p=0 TI=13 TF=28 TF2=29 Meson=h385_s FitTwoPoint
fi

# Simultaneous Kaon fits, multiple momenta
if [ -v DoDisp ]; then
  InDir=$PlotData/corr/2ptp2
  OutDir=$Ensemble/MELFit/2ptp2/s_l

  aFitFiles=()
  for((i=0; i <= MaxPSq; ++i)); do
    aFitFiles+=($InDir/s_l_p2_${i}_g5P_g5P.fold.$DataSeed)
  done
  MultiFit="MultiFit --Hotelling 0 --overwrite --debug-signals --strict 3"
  [ -v FitOptions ] && MultiFit="$MultiFit $FitOptions"

  # Compare unthinned with a couple of choices of thinning
  aTimes=(7:31 10:26 9:21 9:21 9:19 7:28 7:23)
  for (( i=0; i<1; ++i )); do
  (
    case $i in # Not sure ExtraName a good idea? Not used yet anyway
      2) aThin=(1:5:2 1:1:3 1:1:4 4 4 4 4); ExtraName=Test;;
      1) aThin=(1:5:2 1:1:3 1:1:3 1:1:3 1:1:3 1:1:3 1:1:3); ExtraName=Alt;;
      *) aThin=(1:1:2 1:1:3 1:1:3 1:1:3 1:1:3 1:1:3 1:1:3);;
    esac
    yrange=0.187:0.45 SimulP s_l # Simultaneous fits of point-point data at all momenta
    yrange=0.187:0.45 ExtraName=continuum FitOptions='--dispersion continuum' SimulP s_l
  )
  done
fi

# This is a test to see whether the axial at n^2=3 really is compatible with pseudoscalar
# Conclusion: Yes it is - to within statistics.
# We get same range of fitted energies with early/late wall fit on pseudoscalar only
if [ -v DoAxial ]; then
(
  Stat=0
  InPrefix=$PlotData/corr/2ptp2/s_l_p2_3_
  OutPrefix=MELFit/Multi/2ptp2/
  RefText="s-l (n^2=3) E_0="
  ti='5 5 5 5'
  tf='30 30 30 30'
  yrange=0.28:0.335
  NoEnh=

  # Pseudoscalar point pPseudoscalar wall
  FileList=(g5P_g5P g5P_g5W)
  FileArgs=(t=9:21 t=8:21,e=1)
  FitMulti

  # Pseudoscalar point pPseudoscalar wall
  FileList=(g5P_g5P g5P_g5W)
  FileArgs=(t=9:21 t=12:21,e=1)
  FitMulti

  # Axial pt
  FileList=(gT5P_g5P)
  FileArgs=(t=5:17)
  FitMulti

  # Axial pt 1-exp
  FileList=(gT5P_g5P)
  FileArgs=(t=12:21)
  NumExp=1 FitMulti

  # Axial wall 1-exp
  FileList=(gT5P_g5W)
  FileArgs=(t=13:21)
  NumExp=1 FitMulti

  # Axial pt and Axial wall 1-exp
  FileList=(gT5P_g5P gT5P_g5W)
  FileArgs=(t=12:21 t=13:21)
  NumExp=1 FitMulti

  # Axial pt and Axial wall
  FileList=(gT5P_g5P gT5P_g5W)
  FileArgs=(t=5:17 t=13:21,e=1)
  FitMulti

  # Pseudoscalar pt and axial pt
  FileList=(g5P_g5P gT5P_g5P)
  FileArgs=(t=9:21 t=7:21)
  FitMulti

  # All four
  FileList=(g5P_g5P gT5P_g5P g5P_g5W gT5P_g5W)
  FileArgs=(t=9:21 t=7:21 t=9:11,e=1 t=9:11,e=1)
  FitMulti

  # All four
  FileList=(g5P_g5P gT5P_g5P g5P_g5W gT5P_g5W)
  FileArgs=(t=9:21 t=7:21 t=12:21,e=1 t=12:21,e=1)
  FitMulti
)
fi
