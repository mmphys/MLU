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
  x=EL plotglobfit.sh $Dir/{Simul,Separate}/renorm**_all/*.model_fit.txt
  x=qSq xrange=-0.1:1.75 plotglobfit.sh $Dir/{Simul,Separate}/renorm**_all/*_fplus.*.model_fit.txt
  x=qSq xrange=-0.1:2.2 plotglobfit.sh $Dir/{Simul,Separate}/renorm**_all/*_f0.*.model_fit.txt
else
  export MLUf0ELMin=5.1679096878011405468e+08
  export MLUf0ELMax=1.0829813964236147404e+09
  export MLUfplusELMin=6.2578971496614301205e+08
  export MLUfplusELMax=1.0829813964236147404e+09

  Simul= Disable=035 FitSeries=renorm_mostly ContFit.sh&
  Simul= Disable=035 ContFit.sh&
  Simul= Disable=35 ContFit.sh&
  Simul= Disable=35 DisableP=1 ContFit.sh&

  Disable=035 ContFit.sh&
  Disable=35 ContFit.sh&
  Disable=35 DisableP=1 ContFit.sh& # D ) Omit (E/Lambda)^3
  FitOptions=--noblock DisableZ=35 ContFit.sh& # This is my reference fit
  #FitOptions=--noblock DisableZ=2 ContFit.sh&
  #FitOptions=--noblock DisableZ=3 ContFit.sh&
  #FitOptions=--noblock DisableZ=5 ContFit.sh&
  #FitOptions=--noblock DisableZ=23 ContFit.sh&
  #FitOptions=--noblock DisableZ=25 ContFit.sh&

  # These are other variations
  OutDir=Omit FitOptions=--noblock Disable=0 DisableZ=35 ContFit.sh& # A) Omit FV
  OutDir=Omit FitOptions=--noblock Disable=4 DisableZ=35 ContFit.sh& # B) Omit a Lambda
  OutDir=Omit FitOptions=--noblock Disable=1 DisableZ=35 ContFit.sh& # C) Omit Delta Mpi / Lambda
  #OutDir=Omit FitOptions=--noblock Disable=35 DisableP=1 ContFit.sh& # D ) Omit (E/Lambda)^3
  OutDir=Omit FitOptions=--noblock Disable=5 DisableZ=3 ContFit.sh& # D ) Omit (E/Lambda)^3
  EnsOpt=C2 FitOptions=--noblock DisableZ=35 ContFit.sh& # E) Omit Ensemble C2
  EnsOpt=M3 FitOptions=--noblock DisableZ=35 ContFit.sh& # F) Omit Ensemble M3
  EnsOpt=CM FitOptions=--noblock DisableZ=35 ContFit.sh& # G) Omit Ensemble F1M
  EnsOpt=CF FitOptions=--noblock DisableZ=35 ContFit.sh& # H) Omit Ensemble M1, M2, M3
  EnsOpt=MF FitOptions=--noblock DisableZ=35 Disable=4 ContFit.sh& # I) Omit Ensemble C1, C2

  wait
fi
