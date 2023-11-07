#!/usr/bin/env bash

#set -x
set -e

###################################################

# User input

# Optional
#  D:         Number of discretisation constants (default: 1)
#  E:         Number of energy constants (default: 3)
#  Separate:  Set to anything to perform separate fits for f0 and fplus
#  Some:      Set to anything to drop the highest n^2 data points from C ensembles
#  FitSeries: Name of series to fit
#  Series:    List of 'Ensemble seriesname' pairs to override specific choices
#  Disable:   List of constants to disable (each between 0 and 4) both form factors
#  DisableZ:  List of constants to disable (each between 0 and 4) f0 only
#  DisableP:  List of constants to disable (each between 0 and 4) f+ only
#  OutDir:    Where to put the output
#  EnsOpt:    Ensemble options

###################################################

Separate=$((0${Separate+1}))
Some=$((0${Some+1}))
FitSeries=${FitSeries:-renorm}
#OutDir=${OutDir:-C1C2F1MM1M2M3}

eval declare -A ASeries=($Series)

if ! [ -v EnsOpt ]; then
  EnsembleList="C1 C2 F1M M1 M2 M3"
else
  case "$EnsOpt" in
    C1 ) EnsembleList="C2 F1M M1 M2 M3";;
    C2 ) EnsembleList="C1 F1M M1 M2 M3";;
    M1 ) EnsembleList="C1 C2 F1M M2 M3";;
    M2 ) EnsembleList="C1 C2 F1M M1 M3";;
    M3 ) EnsembleList="C1 C2 F1M M1 M2";;
    F1M ) EnsembleList="C1 C2 M1 M2 M3";;
    MF ) EnsembleList="M1 M2 M3 F1M";;
    CF ) EnsembleList="C1 C2 F1M";;
    *) "Unrecognised EnsOpt \"$EnsOpt\""; exit 1;;
  esac
  [ -v OutDir ] && OutDir+=_
  OutDir+=$EnsOpt
fi

if [ -v OutDir ]; then
    ((Separate)) && OutDir+=_Sep
  else
    ((Separate)) && OutDir=Separate || OutDir=Simul
fi

###################################################

# Derived variables

###################################################

unset Disabled
[ -v D ] && Disabled+="D$D"
[ -v E ] && Disabled+="E$E"
if [ "$Disable" != '' ] || [ "$DisableZ" != '' ] || [ "$DisableP" != '' ]
then
  Disabled+="-C$Disable"
  [ "$DisableZ" != '' ] && Disabled+="Z$DisableZ"
  [ "$DisableP" != '' ] && Disabled+="P$DisableP"
  DisableZ="$Disable$DisableZ"
  DisableP="$Disable$DisableP"
fi
OutDir=Cont/$OutDir/${FitSeries}${Disabled}

OutPrefix="F3_K_Ds.corr_"
OutModel=".g5P_g5W.model"

declare -A EnsemblePMax
EnsemblePMax=(C1 4 C2 4 F1M 6 M1 4 M2 4 M3 4)

###################################################

#  Set the list of Files

###################################################

function GetFiles()
{
  unset Files # Returned
  local pMax Ens i ThisSeries
  for Ens in $EnsembleList
  do
    ThisSeries=${ASeries[$Ens]}
    [ -z "$ThisSeries" ] && ThisSeries=$FitSeries
    pMax=${EnsemblePMax[$Ens]}
    if [ "${Ens:0:1}" == C ]; then pMax=$((pMax-Some)); fi
    [ -v Files ] && Files+=' '
    Files+="$Ens/FormFactor/$ThisSeries/3sm_sp2/'*_p2_[0-$pMax].g*'.h5"
  done
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
  local Cmd="$Continuum -o $OutDir/"
  local LogBase="$OutDir/$OutPrefix$FFS$OutModel"
  local ff File FFSwitch Args
  mkdir -p $OutDir
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
[ -v D ] && Continuum+=" -d $D"
[ -v E ] && Continuum+=" -e $E"
Continuum+=" --overwrite"
Continuum+=" --model ${MLUCache}EnsembleInfo.h5"
Continuum+=" -i $HOME/NoSync/"
[ -v FitOptions ] && Continuum+=" $FitOptions"

  GetFiles
  if ((Some)); then
    OutDir+=_some
    qSqRangeMin='*'
  else
    qSqRangeMin='-0.1'
  fi
  if ! [ -v PlotOnly ]; then
  if (( Separate )); then
    DoFit f0
    DoFit fplus
  else
    DoFit f0_fplus
  fi
  fi
  # Now plot it
  for ff in f0 fplus
  do
    [ $ff == fplus ] && qSqRangeMax='1.75' || qSqRangeMax='2.2'
    LogBase="$OutDir/$OutPrefix$ff$OutModel"
    (( Separate )) || KillLogBase
    DoCmd "x=EL xrange='0.48:*' KeyOffset='char 2,0' plotglobfit.sh '$LogBase.txt'"
    DoCmd "x=qSq xrange='$qSqRangeMin:$qSqRangeMax' plotglobfit.sh '$LogBase.txt'"
  done
