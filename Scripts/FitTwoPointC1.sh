#!/usr/bin/env bash

# Simple two-point plots

#set -x
set -e
export Ensemble=C1

if [ "${DoAll+x}" = x ]; then

# D_s
p=0 TI=6 TF=29 NumExp=2 TI2=5 Meson=h6413_s yrange='1.08:1.12' ti='8 8' FitTwoPoint.sh # Preferred
p=0 TI=6 TF=29 NumExp=2 TI2=13 TF2=27 NumExp2=1  Meson=h6413_s yrange='1.08:1.12' ti='8 8' FitTwoPoint.sh
p=0 TI=6 TF=29 NumExp=2 TI2=17 TF2=27 NumExp2=1  Meson=h6413_s yrange='1.08:1.12' ti='8 8' FitTwoPoint.sh

# Kaon
p=0 TI=6 TF=22 NumExp=2 TI2=5 TF2=23 Meson=s_l yrange='0.3:0.315' NumExp2=2 ti='3 3' LabelTF=30 FitTwoPoint.sh
p=0 TI=6 TF=22 NumExp=2 TI2=7 TF2=23 Meson=s_l yrange='0.3:0.315' NumExp2=1 ti='3 3' LabelTF=30 FitTwoPoint.sh
p=0 TI=6 TF=22 NumExp=2 TI2=12 TF2=23 Meson=s_l yrange='0.3:0.315' NumExp2=1 ti='3 3' LabelTF=30 FitTwoPoint.sh


p=1 TI=6 TF=22 NumExp=2 NumExp2=2 TI2=5  TF2=23 Meson=s_l yrange='0.37:0.44' ti='3 3' LabelTF=30 FitTwoPoint.sh
p=1 TI=6 TF=14 NumExp=2 NumExp2=2 TI2=5  TF2=23 Meson=s_l yrange='0.37:0.44' ti='3 3' LabelTF=30 FitTwoPoint.sh
p=1 TI=6 TF=14 NumExp=2 NumExp2=1 TI2=12 TF2=23 Meson=s_l yrange='0.37:0.44' ti='3 3' LabelTF=30 FitTwoPoint.sh # My preferred
p=1 TI=6 TF=22 NumExp=2 NumExp2=1 TI2=12 TF2=23 Meson=s_l yrange='0.37:0.44' ti='3 3' LabelTF=30 FitTwoPoint.sh

for TI2 in 5 12
do
  p=2 TI=6 TF=22 NumExp=2 TF2=23 NumExp2=1 Meson=s_l yrange='0.43:0.5' FitTwoPoint.sh
  p=3 TI=6 TF=22 NumExp=2 TF2=23 NumExp2=1 Meson=s_l yrange='0.48:0.6' FitTwoPoint.sh
  p=4 TI=6 TF=22 NumExp=2 TF2=23 NumExp2=1 Meson=s_l yrange='0.48:0.8' FitTwoPoint.sh
done
fi


#Mtg Thu 19 Jan
p=0 TI=9 TF=29 NumExp=2 TI2=7 TF2=29 NumExp2=2  Meson=h6413_s yrange='1.08:1.12' ti='8 8' FitTwoPoint.sh
#End mtg
