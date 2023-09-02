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
  FitOptions=--block Separate= Disable=035 FitSeries=renorm_mostly ContFit.sh& # Ch 2
  FitOptions=--block Separate= Disable=035 ContFit.sh& # Ch 3
  FitOptions=--block Separate= Disable=35 ContFit.sh& # Ch 4
  FitOptions=--block Separate= Disable=35 DisableP=1 ContFit.sh& # Ch 5
  FitOptions=--block Disable=035 ContFit.sh& # Ch 6
  FitOptions=--block Disable=35 ContFit.sh& # Ch 7
  FitOptions=--block Disable=35 DisableP=1 ContFit.sh& # Ch 8

  #DisableZ=2 ContFit.sh&
  #DisableZ=3 ContFit.sh&
  #DisableZ=5 ContFit.sh&
  #DisableZ=23 ContFit.sh&
  #DisableZ=25 ContFit.sh&
  DisableZ=35 ContFit.sh& # This is my reference fit

  # These are other variations
  OutDir=Omit Disable=0 DisableZ=35 ContFit.sh& # A) Omit FV
  OutDir=Omit Disable=4 DisableZ=35 ContFit.sh& # B) Omit a Lambda
  OutDir=Omit Disable=1 DisableZ=35 ContFit.sh& # C) Omit Delta Mpi / Lambda
  OutDir=Omit Disable=5 DisableZ=3 ContFit.sh& # D ) Omit (E/Lambda)^3
  EnsOpt=C2 DisableZ=35 ContFit.sh& # E) Omit Ensemble C2
  EnsOpt=M3 DisableZ=35 ContFit.sh& # F) Omit Ensemble M3
  EnsOpt=CM DisableZ=35 ContFit.sh& # G) Omit Ensemble F1M
  EnsOpt=CF DisableZ=35 ContFit.sh& # H) Omit Ensemble M1, M2, M3
  EnsOpt=MF DisableZ=35 Disable=4 ContFit.sh& # I) Omit Ensemble C1, C2

  wait
fi
