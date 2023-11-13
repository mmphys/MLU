#!/usr/bin/env bash

#set -x
set -e

###################################################

# Perform all the continuum fit variations

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
else
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

  # Historic fits documented in Continuum.tex
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
  EnsOpt=C1 DisableZ=34 ContFit.sh& # b) Omit Ensemble C1
  EnsOpt=C2 DisableZ=34 ContFit.sh& # c) Omit Ensemble C2
  EnsOpt=M1 DisableZ=34 ContFit.sh& # d) Omit Ensemble M1
  EnsOpt=M2 DisableZ=34 ContFit.sh& # e) Omit Ensemble M2
  EnsOpt=M3 DisableZ=34 ContFit.sh& # f) Omit Ensemble M3
  EnsOpt=CF DisableZ=34 ContFit.sh& # g) Omit Ensemble M1, M2, M3
  Some= DisableZ=34 ContFit.sh& # h) Omit n^2_max from C1 C2
  OutDir=Pole300-25 DisableZ=34 FitOptions='--poles=3e8 --polev=-2.5e7' ContFit.sh& # i) Manual pole
  NameExtra=PoleV DisableZ=34 FitOptions='--polev=-50e6' ContFit.sh& # j) Manual pole
  OutDir=PoleS250 DisableZ=34 FitOptions='--poles=25e7' ContFit.sh& # k) Manual pole
  OutDir=AltC1 Series='C1 renormold' DisableZ=34 ContFit.sh& # l) Alternate C1 fits

  # These are destructive tests - not part of my error budget
  EnsOpt=F1M DisableZ=34 ContFit.sh& # m) Omit Ensemble F1M
  EnsOpt=MF DisableZ=34 D=0 ContFit.sh& # n) Omit Ensemble C1, C2
  OutDir=Omit Disable=X DisableZ=34 ContFit.sh& # o) Omit Chiral
  OutDir=Omit Disable=VX E=2 DisableZ=3 ContFit.sh& # p) Omit FV and Chiral
  OutDir=Omit D=0       DisableZ=34 ContFit.sh& # q) Omit a Lambda
  OutDir=Omit Disable=1 DisableZ=34 ContFit.sh& # r) Omit Delta Mpi / Lambda
  OutDir=Omit E=2 DisableZ=3 ContFit.sh& # s) Omit (E/Lambda)^3
  wait
fi
