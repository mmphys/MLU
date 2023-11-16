#!/usr/bin/env bash

#set -x
set -e

###################################################

# Perform all the continuum fit variations

# Input:
# Plot: set to anything to perform plots only
# Do:   original  - fits documented in Continuum.tex (*Default*)
#       shrink    - newer fits using shrinkage
#       linear    - drop as few data points as possible to make linear fit work without shrinkage
#                 - anything else do nothing
# DataRange:  set to anything so that plots cover range with data only
#             otherwise complete physical range covered

###################################################

Do=${Do-original shrink linear}

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

# Historic fits documented in Continuum.tex
# and the first version in my thesis

###################################################

function SensitivityOriginal()
{
  local AltPoleS='--poles=313e6'
  local AltPoleV='--polev=-5e6'
  FitOptions=--block Separate= E=1 Disable=V FitSeries=renorm_mostly ContFit.sh& # Ch 2
  FitOptions=--block Separate= E=1 Disable=V ContFit.sh& # Ch 3
  FitOptions=--block Separate= E=1 ContFit.sh& # Ch 4
  FitOptions=--block Separate= E=1 DisableP=1 ContFit.sh& # Ch 5
  FitOptions=--block E=1 Disable=V ContFit.sh& # Ch 6
  FitOptions=--block E=1 ContFit.sh& # Ch 7
  FitOptions=--block E=1 DisableP=1 ContFit.sh& # Ch 8

  #DisableZ=2 ContFit.sh&
  #DisableZ=3 ContFit.sh&
  #DisableZ=4 ContFit.sh&
  #DisableZ=23 ContFit.sh&
  #DisableZ=24 ContFit.sh&
  DisableZ=34 ContFit.sh& # Ch 9 This is my reference fit

  # These are other variations
  OutDir=Omit Disable=V DisableZ=34 ContFit.sh& # a) Omit FV
  OutDir=Omit EnsOpt=C1 DisableZ=34 ContFit.sh& # b) Omit Ensemble C1
  OutDir=Omit EnsOpt=C2 DisableZ=34 ContFit.sh& # c) Omit Ensemble C2
  OutDir=Omit EnsOpt=M1 DisableZ=34 ContFit.sh& # d) Omit Ensemble M1
  OutDir=Omit EnsOpt=M2 DisableZ=34 ContFit.sh& # e) Omit Ensemble M2
  OutDir=Omit EnsOpt=M3 DisableZ=34 ContFit.sh& # f) Omit Ensemble M3
  OutDir=Omit EnsOpt=M DisableZ=34 ContFit.sh& # g) Omit Ensemble M1, M2, M3
  Some= DisableZ=34 ContFit.sh& # h) Omit n^2_max from C1 C2
  NameExtra=PoleSV DisableZ=34 FitOptions="$AltPoleS $AltPoleV" ContFit.sh& # i) move S&V poles
  NameExtra=PoleV  DisableZ=34 FitOptions="$AltPoleV" ContFit.sh& # j) Move vector pole
  NameExtra=PoleS  DisableZ=34 FitOptions="$AltPoleS" ContFit.sh& # k) Move scalar pole
  NameExtra=AltC1 Series='C1 renormold' DisableZ=34 ContFit.sh& # Ch 10 l) Alternate C1 fits

  # These are destructive tests - not part of my error budget
  OutDir=Omit EnsOpt=F1M DisableZ=34 ContFit.sh& # m) Omit Ensemble F1M
  OutDir=Omit EnsOpt=C DisableZ=34 D=0 ContFit.sh& # n) Omit Ensemble C1, C2
  OutDir=Omit Disable=X DisableZ=34 ContFit.sh& # o) Omit Chiral
  OutDir=Omit Disable=VX E=2 DisableZ=3 ContFit.sh& # p) Omit FV and Chiral
  OutDir=Omit D=0       DisableZ=34 ContFit.sh& # q) Omit a Lambda
  OutDir=Omit Disable=1 DisableZ=34 ContFit.sh& # r) Omit Delta Mpi / Lambda
  OutDir=Omit E=2 DisableZ=3 ContFit.sh& # s) Omit (E/Lambda)^3
}

###################################################

# New reference fit - linear with shrinkage

###################################################

function SensitivityShrinkLinear()
{
  local AltPoleS='--poles=313e6'
  local AltPoleV='--polev=-5e6'
  local E=1
  export E
  ContFit.sh& # Ref) This is my reference fit
  OutDir=Omit Disable=V ContFit.sh& # a) Omit FV
  OutDir=Omit EnsOpt=C1 ContFit.sh& # b) Omit Ensemble C1
  OutDir=Omit EnsOpt=C2 ContFit.sh& # c) Omit Ensemble C2
  OutDir=Omit EnsOpt=M1 ContFit.sh& # d) Omit Ensemble M1
  OutDir=Omit EnsOpt=M2 ContFit.sh& # e) Omit Ensemble M2
  OutDir=Omit EnsOpt=M3 ContFit.sh& # f) Omit Ensemble M3
  OutDir=Omit EnsOpt=M ContFit.sh& # g) Omit Ensemble M1, M2, M3
  Some= ContFit.sh& # h) Omit n^2_max from C1 C2
  NameExtra=PoleSV FitOptions="$AltPoleS $AltPoleV" ContFit.sh& # i) move S&V poles
  NameExtra=PoleV  FitOptions="$AltPoleV" ContFit.sh& # j) Move vector pole
  NameExtra=PoleS  FitOptions="$AltPoleS" ContFit.sh& # k) Move scalar pole
  NameExtra=AltC1 Series='C1 renormold' ContFit.sh& # Ch 10 l) Alternate C1 fits
  unset E
  unset Shrink
  DisableZ=34 ContFit.sh& # m) Original reference fit
}

# Linear model with shrinkage

function SensitivityShrink()
{
  local Shrink=0.005
  export Shrink
  SensitivityShrinkLinear
}

# Drop as few data points as possible to make linear model work without shrinkage

function SensitivityLinear()
{
  local PMaxFZero='C1 4 C2 3 F1M 6 M1 3 M2 3 M3 3'
  export PMaxFZero
  SensitivityShrinkLinear
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
      original) SensitivityOriginal;;
      shrink)   SensitivityShrink;;
      linear)   SensitivityLinear;;
      *)        echo "Unrecognised fit variation - doing nothing";;
    esac
  done
  wait || : # Because waiting if there are no background jobs will return 127
  )
fi
