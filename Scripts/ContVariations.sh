#!/usr/bin/env bash

#set -x
set -e

###################################################

# Perform all the continuum fit variations

###################################################

if [ -v Plot ]
then
  shopt -s globstar
  Dir=/Users/mike/NoSync/Cont
  x=EL plotglobfit.sh $Dir/**/renorm**/*.model_fit.txt
  x=qSq xrange=-0.1:1.75 plotglobfit.sh $Dir/**/renorm**/*_fplus.*.model_fit.txt
  x=qSq xrange=-0.1:2.2 plotglobfit.sh $Dir/**/renorm**/*_f0.*.model_fit.txt
else
  export MLUf0ELMin=5.1679096878011405468e+08
  export MLUf0ELMax=1.0829813964236147404e+09
  export MLUfplusELMin=6.2578971496614301205e+08
  export MLUfplusELMax=1.0829813964236147404e+09

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
  DisableZ=34 ContFit.sh& # Ch 9 This is my reference fit. # L) _some=Omit n^2_max from C1 C2

  # These are other variations
  OutDir=Omit Disable=V DisableZ=34 ContFit.sh& # A) Omit FV
  OutDir=Omit Disable=X DisableZ=34 ContFit.sh& # B) Omit Chiral
  OutDir=Omit Disable=VX E=2 DisableZ=3 ContFit.sh& # C) Omit FV and Chiral
  OutDir=Omit D=0       DisableZ=34 ContFit.sh& # D) Omit a Lambda
  OutDir=Omit Disable=1 DisableZ=34 ContFit.sh& # E) Omit Delta Mpi / Lambda
  OutDir=Omit E=2 DisableZ=3 ContFit.sh& # F ) Omit (E/Lambda)^3
  EnsOpt=C2 DisableZ=34 ContFit.sh& # G) Omit Ensemble C2
  EnsOpt=M3 DisableZ=34 ContFit.sh& # H) Omit Ensemble M3
  EnsOpt=CM DisableZ=34 ContFit.sh& # I) Omit Ensemble F1M
  EnsOpt=CF DisableZ=34 ContFit.sh& # J) Omit Ensemble M1, M2, M3
  EnsOpt=MF DisableZ=34 D=0 ContFit.sh& # K) Omit Ensemble C1, C2

  wait
fi
