#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=F1M
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit
fi
. PlotCommon.sh
. FitTwoPoint.sh

declare -A ayRange
ayRange[h385_s,0]='0.72:0.735'
ayRange[s_l,0]='0.187:0.197'
ayRange[s_l,1]='0.225:0.245'
ayRange[s_l,2]='0.25:0.28'
ayRange[s_l,3]='0.28:0.32'
ayRange[s_l,4]='0.3:0.37'
ayRange[s_l,5]='0.32:0.38'
ayRange[s_l,6]='0.2:0.42'

if [ "${DoAll+x}" = x ]; then

# D_s
p=0 TI=6 TF=29 NumExp=3 TI2=5 Meson=h385_s yrange='0.72:0.735' FitTwoPoint

# Kaon
p=0 TI=10 TF=26 NumExp=2 TI2=7 Meson=s_l FitTwoPoint
p=1 TI=8 TF=28 NumExp=2 TI2=9 NumExp2=1 Meson=s_l FitTwoPoint
p=2 TI=8 TF=26 NumExp=2 TI2=8 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
p=3 TI=8 TF=28 NumExp=2 TI2=8 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
p=3 TI=8 TF=28 NumExp=2 TI2=6 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
#p=3 TI=8 TF=28 NumExp=2 TI2=6 TF2=22 NumExp2=1 Meson=s_l FitTwoPoint #also good
p=4 TI=10 TF=28 NumExp=2 TI2=8 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
#p=4 TI=10 TF=25 NumExp=2 TI2=8 TF2=23 NumExp2=1 Meson=s_l FitTwoPoint #also good
p=5 TI=10 TF=28 NumExp=2 TI2=12 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
p=6 TI=10 TF=28 NumExp=2 TI2=12 TF2=28 NumExp2=1 Meson=s_l FitTwoPoint
fi
