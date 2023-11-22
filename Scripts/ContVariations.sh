#!/usr/bin/env bash

#set -x
set -e

###################################################

# Perform all the continuum fit variations

# Input:
# Plot: set to anything to perform plots only
# Do:   cubic     - fits documented in Continuum.tex (*Default*)
#       shrink    - newer fits using shrinkage
#       linear    - drop as few data points as possible to make linear fit work without shrinkage
#                 - anything else do nothing
# DataRange:  set to anything so that plots cover range with data only
#             otherwise complete physical range covered

###################################################

Do=${Do-cubic shrink linear}

# Set domain - only used when performing fits
# Can be useful to source this script when manually doing global fits

if [ -v DataRange ]; then
  # Range covered by the data points
  export MLUf0ELMin=5.1679097012296748161e+08
  export MLUf0ELMax=1.0829813959019107819e+09
  export MLUfplusELMin=6.2578971607197594643e+08
else
  # Physical range
  export MLUf0ELMin=4.95644e+08
  export MLUf0ELMax=1.0829813959019107819e+09
  export MLUfplusELMin=$MLUf0ELMin
fi
export MLUfplusELMax=$MLUf0ELMax

###################################################

# Sensitivity analysis common across different reference fits

###################################################

function SensitivityCommon()
{
  local AltPoleS='--poles=313e6'
  local AltPoleV='--polev=-5e6'
  ContFit.sh& # Ref) This is my reference fit
  OutDir=Omit Disable=V ContFit.sh& # a) Omit FV
  OutDir=Omit EnsOpt=C1 ContFit.sh& # b) Omit Ensemble C1
  OutDir=Omit EnsOpt=C2 ContFit.sh& # c) Omit Ensemble C2
  OutDir=Omit EnsOpt=M1 ContFit.sh& # d) Omit Ensemble M1
  OutDir=Omit EnsOpt=M2 ContFit.sh& # e) Omit Ensemble M2
  OutDir=Omit EnsOpt=M3 ContFit.sh& # f) Omit Ensemble M3
  OutDir=Omit EnsOpt=M ContFit.sh& # g) Omit Ensemble M1, M2, M3
  #Some= ContFit.sh& # h) Omit n^2_max from C1 C2
  PMaxFZero='C1 3 C2 3 F1M 6 M1 3 M2 3 M3 3' \
  PMaxFPlus='C1 3 C2 3 F1M 6 M1 3 M2 3 M3 3' \
  OutDir=Omit NameExtra=NMaxCM ContFit.sh& # h) Omit n^2_max from C & M ensembles
  NameExtra=PoleSV FitOptions="$AltPoleS $AltPoleV" ContFit.sh& # i) move S&V poles
  NameExtra=PoleV  FitOptions="$AltPoleV" ContFit.sh& # j) Move vector pole
  NameExtra=PoleS  FitOptions="$AltPoleS" ContFit.sh& # k) Move scalar pole
  NameExtra=AltC1 Series='C1 renormold' ContFit.sh& # Ch 10 l) Alternate C1 fits
}

###################################################

# Historic fits documented in Continuum.tex
# and the first version in my thesis

###################################################

function SensitivityOriginal()
{
  FitOptions=--block Separate= Disable=V FitSeries=renorm_mostly ContFit.sh& # Ch 2
  FitOptions=--block Separate= Disable=V ContFit.sh& # Ch 3
  FitOptions=--block Separate= ContFit.sh& # Ch 4
  FitOptions=--block Separate= DisableP=1 ContFit.sh& # Ch 5
  FitOptions=--block OutDir=Block Disable=V ContFit.sh& # Ch 6
  FitOptions=--block OutDir=Block ContFit.sh& # Ch 7
  FitOptions=--block OutDir=Block DisableP=1 ContFit.sh& # Ch 8
  ContFit.sh& # Full covariance matrix with linear model - terrible Hotelling p-Value

  local E=3
  export E

  #DisableZ=2 ContFit.sh&
  #DisableZ=3 ContFit.sh&
  #DisableZ=4 ContFit.sh&
  #DisableZ=23 ContFit.sh&
  #DisableZ=24 ContFit.sh&

  local DisableZ=34
  export DisableZ
  SensitivityCommon

  # These are destructive tests - not part of my error budget
  local OutDir=Omit
  export OutDir
  EnsOpt=F1M ContFit.sh& # m) Omit Ensemble F1M
  EnsOpt=C D=0 ContFit.sh& # n) Omit Ensemble C1, C2
  Disable=X ContFit.sh& # o) Omit Chiral
  Disable=VX E=2 DisableZ=3 ContFit.sh& # p) Omit FV and Chiral
  D=0       ContFit.sh& # q) Omit a Lambda
  Disable=1 ContFit.sh& # r) Omit Delta Mpi / Lambda
  E=2 DisableZ=3 ContFit.sh& # s) Omit (E/Lambda)^3
}

###################################################

# New reference fit - linear with shrinkage

###################################################

# Linear model with shrinkage

function SensitivityShrink()
{
  E=3 DisableZ=34 ContFit.sh& # m) Original reference fit
  local Shrink=0.005
  export Shrink
  SensitivityCommon
}

# Drop as few data points as possible to make linear model work without shrinkage

function SensitivityLinear()
{
  E=3 DisableZ=34 ContFit.sh& # m) Original reference fit
  local PMaxFZero='C1 4 C2 3 F1M 6 M1 3 M2 3 M3 3'
  export PMaxFZero
  SensitivityCommon
}

###################################################

# Main loop

###################################################

if [ -v Plot ]
then
  (
  shopt -s globstar nullglob
  Dir=/Users/mike/NoSync/Cont
  for f in $Dir/**/renorm**/*.model_fit.txt; do
    x=EL xrange='0.48:*' KeyOffset='char 2,0' plotglobfit.sh "$f" &
  done
  for f in $Dir/**/renorm**/*_fplus.*.model_fit.txt; do
    x=qSq xrange=-0.1:1.75 plotglobfit.sh "$f" &
  done
  for f in $Dir/**/renorm**/*_f0.*.model_fit.txt; do
    x=qSq xrange=-0.1:2.2 plotglobfit.sh "$f" &
  done
  wait
  )
elif [ -n "$Do" ]
then
  echo "Making fit variations $Do"
  (
  export Do
  for Do in ${Do,,}
  do
    case ${Do,,} in
      cubic)    SensitivityOriginal;;
      shrink)   SensitivityShrink;;
      linear)   SensitivityLinear;;
      *)        echo "Unrecognised fit variation - doing nothing";;
    esac
  done
  wait || : # Because waiting if there are no background jobs will return 127
  )
fi
