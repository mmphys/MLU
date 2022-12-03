#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=F1M

# D_s
p=0 TI=6 TF=29 NumExp=3 TI2=5 Meson=h385_s yrange='0.72:0.735' FitTwoPoint.sh

# These are the best two Kaon fits, working with Peter Friday 25 Nov 2022
#p=0 TI=10 TF=26 NumExp=2 TI2=7 Meson=s_l yrange='0.187:0.197' FitTwoPoint.sh
#p=1 TI=8 TF=28 NumExp=2 TI2=9 NumExp2=1 Meson=s_l yrange='0.225:0.245' FitTwoPoint.sh
p=2 TI=8 TF=26 NumExp=2 TI2=8 TF2=28 NumExp2=1 Meson=s_l yrange='0.250:0.280' FitTwoPoint.sh
