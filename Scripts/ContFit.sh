#!/usr/bin/env bash

#set -x
set -e

# User input
# Disable: space separated list of constants to disable (between 0 and 4)

FitSeries=${FitSeries:-renorm}

Enable=(1 1 1 1 1)
unset ModelArgs
unset Disabled
for MyDisable in $Disable; do
  Enable[$MyDisable]=0
  ModelArgs=${ModelArgs-,EnableC${MyDisable}=false}
  Disabled+=$MyDisable
done
OutDir=Cont/${OutDir:-C1C2F1MM1M2M3}/${FitSeries}${Disabled:+-C${Disabled}}_

# Now perform fit

Continuum="Continuum --debug-signals --Hotelling 0 --summary 2"
Continuum+=" --overwrite"
Continuum+=" --model ${MLUCache}EnsembleInfo.h5"
Continuum+=" -i $HOME/NoSync/"

pMax=4
pMaxF1M=6
xRangeMin='-0.1'
for Some in 0 1
do
  pString="3sm_sp2/'*_p2_[0-$((pMax-Some))].g*'.h5"
  pStringF1M="3sm_sp2/'*_p2_[0-$((pMaxF1M-Some))].g*'.h5"
  Files="C1/FormFactor/$FitSeries/$pString$ModelArgs"
  Files="$Files C2/FormFactor/$FitSeries/$pString$ModelArgs"
  Files="$Files F1M/FormFactor/renorm/$pStringF1M$ModelArgs"
  Files="$Files M1/FormFactor/$FitSeries/$pString$ModelArgs"
  Files="$Files M2/FormFactor/$FitSeries/$pString$ModelArgs"
  Files="$Files M3/FormFactor/$FitSeries/$pString$ModelArgs"
  ((Some)) && OutSubDir="${OutDir}some" || OutSubDir="${OutDir}all"
  mkdir -p $OutSubDir
  LogPrefix="${OutSubDir}/F3_K_Ds.corr_"
  for fplus in 0 1
  do
    Cmd="$Continuum -o $OutSubDir/"
    if ((fplus)); then
      Cmd="$Cmd -f fplus"
      LogBase="${LogPrefix}fplus"
      xRangeMax='1.75'
    else
      LogBase="${LogPrefix}f0"
      xRangeMax='2.2'
    fi
    LogBase="$LogBase.g5P_g5W.model"
    Cmd="$Cmd $Files"
    #echo $Cmd
    echo $Cmd &>  $LogBase.log
    if ! eval $Cmd &>> $LogBase.log; then echo "Returned ${PIPESTATUS[0]}" &>> $LogBase.log; fi
    # Now plot it
    export PlotField x
    for x in '' EL; do
    for PlotField in '' adjusted; do
    case $x in
      EL) unset Cmd;;
      *)  Cmd="xrange='$xRangeMin:$xRangeMax' ";;
    esac
    Cmd+="plotglobfit.sh '$LogBase'"
    echo $Cmd &>> $LogBase.log
    eval $Cmd &>> $LogBase.log
    done
    done
  done
  xRangeMin='*'
done
