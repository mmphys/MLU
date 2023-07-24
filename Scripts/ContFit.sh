#!/usr/bin/env bash

#set -x
set -e

Renorm=${Renorm:-C1C2F1MM1M2M3}
OutDir=Cont/$Renorm
FitSeries=disp #${series:-disp}
MLUSeed=2263212701

Continuum="Continuum --debug-signals --Hotelling 0 --summary 2"
Continuum+=" --overwrite"
Continuum+=" --model ${MLUCache}EnsembleInfo.h5"
Continuum+=" -i $HOME/NoSync/"

pMax=4
pMaxF1M=6
xrange='-0.1:2.2'
for Some in 0 1
do
  pString="3sm_sp2/'*_p2_[0-$((pMax-Some))].g*'.h5"
  pStringF1M="3sm_sp2/'*_p2_[0-$((pMaxF1M-Some))].g*'.h5"
  Files="C1/FormFactor/$FitSeries/$pString"
  Files="$Files C2/FormFactor/$FitSeries/$pString"
  Files="$Files F1M/FormFactor/renorm/$pStringF1M"
  Files="$Files M1/FormFactor/$FitSeries/$pString"
  Files="$Files M2/FormFactor/$FitSeries/$pString"
  Files="$Files M3/FormFactor/$FitSeries/$pString"
  OutSubDir="$OutDir/${FitSeries}_"
  ((Some)) && OutSubDir="${OutSubDir}some" || OutSubDir="${OutSubDir}all"
  mkdir -p $OutSubDir
  LogPrefix="${OutSubDir}/F3_K_Ds.corr_"
  for fplus in 0 1
  do
    Cmd="$Continuum -o $OutSubDir/"
    if ((fplus)); then
      Cmd="$Cmd -f fplus"
      LogBase="${LogPrefix}fplus"
    else
      LogBase="${LogPrefix}f0"
    fi
    LogBase="$LogBase.g5P_g5W.model"
    Cmd="$Cmd $Files"
    #echo $Cmd
    echo $Cmd &>  $LogBase.log
    if ! eval $Cmd &>> $LogBase.log; then echo "Returned ${PIPESTATUS[0]}" &>> $LogBase.log; fi
    # Now plot it
    Cmd="xrange='$xrange' plotglobfit.sh '$LogBase.$MLUSeed.dat'"
    echo $Cmd &>> $LogBase.log
    eval $Cmd &>> $LogBase.log
  done
  xrange='*:2.2'
done
