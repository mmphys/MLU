#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=F1M

if [ "${DoAll+x}" = x ]; then

# D_s
p=0 TI=6 TF=29 NumExp=3 TI2=5 Meson=h385_s yrange='0.72:0.735' FitTwoPoint.sh

# Kaon
p=0 TI=10 TF=26 NumExp=2 TI2=7 Meson=s_l yrange='0.187:0.197' FitTwoPoint.sh
p=1 TI=8 TF=28 NumExp=2 TI2=9 NumExp2=1 Meson=s_l yrange='0.225:0.245' FitTwoPoint.sh
p=2 TI=8 TF=26 NumExp=2 TI2=8 TF2=28 NumExp2=1 Meson=s_l yrange='0.250:0.280' FitTwoPoint.sh
p=3 TI=8 TF=28 NumExp=2 TI2=6 TF2=28 NumExp2=1 Meson=s_l yrange='0.28:0.32' FitTwoPoint.sh
#p=3 TI=8 TF=28 NumExp=2 TI2=6 TF2=22 NumExp2=1 Meson=s_l yrange='0.28:0.32' FitTwoPoint.sh #also good
p=4 TI=10 TF=28 NumExp=2 TI2=8 TF2=28 NumExp2=1 Meson=s_l yrange='0.3:0.37' FitTwoPoint.sh
#p=4 TI=10 TF=25 NumExp=2 TI2=8 TF2=23 NumExp2=1 Meson=s_l yrange='0.3:0.37' FitTwoPoint.sh #also good
p=5 TI=10 TF=28 NumExp=2 TI2=12 TF2=28 NumExp2=1 Meson=s_l yrange='0.32:0.38' FitTwoPoint.sh
p=6 TI=10 TF=28 NumExp=2 TI2=12 TF2=28 NumExp2=1 Meson=s_l yrange='0.2:0.42' FitTwoPoint.sh
fi
