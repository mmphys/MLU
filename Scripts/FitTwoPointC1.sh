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
  p=3 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:23:10:1,e=2 t=7:23,e=1
  p=4 Meson=s_l SubDir=e2_1a TwoPointScan.sh t=3:23:10:1,e=2 t=7:23,e=1
fi

if [ -v DoScanW ]; then
  Meson=h6413_s SubDir=e2_1b TwoPointScan.sh t=6:27,e=2 t=13:27:10:1,e=1
  p=0 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=6:23,e=2 t=3:23:13:1,e=1
  p=1 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=6:23,e=2 t=3:23:13:1,e=1
  p=2 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=5:20,e=2 t=3:20:13:1,e=1
  p=3 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=6:23,e=2 t=3:23:13:1,e=1
  p=4 Meson=s_l SubDir=e2_1b TwoPointScan.sh t=6:23,e=2 t=3:23:13:1,e=1
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
  p=2 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 Meson=s_l yrange='0.42:0.52' FitTwoPoint.sh
  p=3 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 Meson=s_l yrange='0.48:0.6' FitTwoPoint.sh
  p=4 TI=6 TF=23 NumExp=2 TI2=7 TF2=23 NumExp2=1 Meson=s_l yrange='0.48:0.8' FitTwoPoint.sh
fi

