#!/usr/bin/env bash

# Simple two-point plots

set -e

############################################################

# Inputs

############################################################

#Ensemble=${Ensemble:-F1M}
. PlotCommon.sh
PCMakeScriptOutDir
#set -x

OptP="$1" # t=...[,e=...]
OptW="$2"

qHeavy=${qHeavy:-h$Heavy}
Meson=${Meson:-${qHeavy}_s}
p=${p:-0}
FitOptions="--debug-signals --mindp 2 --iter 100000 --strict 3 $FitOptions"

# SubDir - which subdirectory to place output in
# SummaryOptions - extra options for the Summary

MLUSeed=${MLUSeed:-1835672416}
Corr=${Corr:-corr}

############################################################

# Error check

############################################################

unset DoErrorMsg
if [ -z "$OptP" ]
then
  echo "Point options unspecified. Expect as a minimum t=ti:tf"
  DoErrorMsg=
fi
if [ -z "$OptW" ]
then
  echo "Wall options unspecified. Expect as a minimum t=ti:tf"
  DoErrorMsg=
fi
if [ -v DoErrorMsg ]; then exit 1; fi

############################################################

# Derived

############################################################

if [[ ! $p =~ _ ]]; then p="p2_$p"; fi # Assume momentum is p^2 (but allow it to be set explicitly)

InDir=$PlotData/$Corr/2ptp2
OutDir=$Meson
if [ -v SubDir ]; then OutDir="$OutDir/$SubDir"; fi

Suffix=fold.$MLUSeed.h5
FileBase=${Meson}_${p}_g5P_g5
LogFile=$OutDir/${FileBase}W.log

############################################################

# Perform the fit

############################################################

echo "$Ensemble: Scan two point ${Meson}_${p}, Point $OptP, Wall $OptW"

mkdir -p $OutDir # So I can redirect output here

Cmd="MultiFit -o \"$OutDir/\" -i \"$InDir/$FileBase\" $FitOptions"
Cmd="$Cmd P.$Suffix,$OptP"
Cmd="$Cmd W.$Suffix,$OptW"
echo "$Cmd" >> $LogFile
eval  $Cmd &>> $LogFile

Cmd="FitSummary --all --strict 0 -i \"$OutDir/\" -o \"$OutDir/\"" # --errdig 2 --tablen
[ -v SummaryOptions ] && Cmd="$Cmd $SummaryOptions"
Cmd="$Cmd '${Meson}_${p}.*.h5'"
echo "$Cmd" >> $LogFile
eval  $Cmd &>> $LogFile

Cmd="rm \"$OutDir/${Meson}_${p}.corr_\"*"
echo "$Cmd" >> $LogFile
eval  $Cmd &>> $LogFile
