#!/usr/bin/env bash

#set -x
set -e

###################################################

# User input

# Optional
#  Simul:     Set to anything to disable simultaneous fit of f0 and fplus
#  FitSeries: Name of series to fit
#  Disable:   List of constants to disable (each between 0 and 4) both form factors
#  DisableZ:  List of constants to disable (each between 0 and 4) f0 only
#  DisableP:  List of constants to disable (each between 0 and 4) f+ only
#  OutDir:    Where to put the output

###################################################

Simul=$((1-0${Simul+1}))
FitSeries=${FitSeries:-renorm}
#OutDir=${OutDir:-C1C2F1MM1M2M3}
if ! [ -v OutDir ]; then
  ((Simul)) && OutDir=Simul || OutDir=Separate
fi

###################################################

# Derived variables

###################################################

unset Disabled
if [ "$Disable" != '' ] || [ "$DisableZ" != '' ] || [ "$DisableP" != '' ]
then
  Disabled="-C$Disable"
  [ "$DisableZ" != '' ] && Disabled+="Z$DisableZ"
  [ "$DisableP" != '' ] && Disabled+="P$DisableP"
  DisableZ="$Disable$DisableZ"
  DisableP="$Disable$DisableP"
fi
OutDir=Cont/$OutDir/${FitSeries}${Disabled}_

OutPrefix="F3_K_Ds.corr_"
OutModel=".g5P_g5W.model"

###################################################

#  Set the list of Files

###################################################

function GetFiles()
{
  local Some=$1
  local pMax=4
  local pMaxF1M=6
  local pString="3sm_sp2/'*_p2_[0-$((pMax-Some))].g*'.h5"
  local pStringF1M="3sm_sp2/'*_p2_[0-$((pMaxF1M-Some))].g*'.h5"
  Files="C1/FormFactor/$FitSeries/$pString"
  Files+=" C2/FormFactor/$FitSeries/$pString"
  Files+=" F1M/FormFactor/$FitSeries/$pStringF1M"
  Files+=" M1/FormFactor/$FitSeries/$pString"
  Files+=" M2/FormFactor/$FitSeries/$pString"
  Files+=" M3/FormFactor/$FitSeries/$pString"
}

###################################################

# Execute Cmd

###################################################

function KillLogBase()
{
  if [ -e $LogBase.log ]
  then
    rm $LogBase.log
  fi
}

function DoCmd()
{
  local Cmd="$1"
      #echo $Cmd
       echo $Cmd &>> $LogBase.log
  if ! eval $Cmd &>> $LogBase.log
  then
    echo "Returned ${PIPESTATUS[0]}" &>> $LogBase.log
  fi
}

###################################################

# Do one fit

###################################################

function DoFit()
{
  local FFS="$1"
  local Cmd="$Continuum -o $OutSubDir/"
  local LogBase="$OutSubDir/$OutPrefix$FFS$OutModel"
  local ff File FFSwitch Args
  mkdir -p $OutSubDir
  KillLogBase
  for ff in ${FFS//_/ }; do
    FFSwitch+=${FFSwitch:+,}$ff,
    [ "$ff" = "f0" ] && FFSwitch+=$DisableZ || FFSwitch+=$DisableP
    for File in $Files; do
      Args+="${Args:+ }$File,ff=$ff"
    done
  done
  DoCmd "$Cmd -f $FFSwitch $Args"
}

###################################################

# Now perform fit

###################################################

Continuum="Continuum --debug-signals --Hotelling 0 --summary 2"
Continuum+=" --overwrite"
Continuum+=" --model ${MLUCache}EnsembleInfo.h5"
Continuum+=" -i $HOME/NoSync/"

for Some in 0 1
do
  GetFiles $Some
  ((Some)) && OutSubDir="${OutDir}some" || OutSubDir="${OutDir}all"
  ((Some)) && qSqRangeMin='*' || qSqRangeMin='-0.1'
  if (( Simul )); then
    DoFit f0_fplus
  else
    DoFit f0
    DoFit fplus
  fi
  # Now plot it
  for ff in f0 fplus
  do
    [ $ff == fplus ] && qSqRangeMax='1.75' || qSqRangeMax='2.2'
    LogBase="$OutSubDir/$OutPrefix$ff$OutModel"
    (( Simul )) && KillLogBase
    DoCmd "x=EL plotglobfit.sh '$LogBase.txt'"
    DoCmd "x=qSq xrange='$qSqRangeMin:$qSqRangeMax' plotglobfit.sh '$LogBase.txt'"
  done
done