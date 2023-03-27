#!/usr/bin/env bash

# Functions for different types of two-point fits (source me)

############################################################

# Execute command

# Input
#  Cmd:     Command to execute
#  LogFile: Where to redirect output

############################################################

function DoCmd()
{
  #echo $Cmd
  echo $Cmd &>> $LogFile
  if ! eval $Cmd &>> $LogFile; then
    echo "Error ${PIPESTATUS[0]} executing $Cmd"
  fi
}

############################################################

# Perform a simple, two-point fit of point-point and point-wall

############################################################

function FitTwoPoint()
{
########################################
# Input parameters
local qHeavy=${qHeavy:-h$Heavy}
local Meson=${Meson:-${qHeavy}_s}
local p=${p:-0}
local NumExp=${NumExp:-3}
local NumExp2=${NumExp2:-${NumExp}}
local MaxExp=$(( NumExp2 > NumExp ? NumExp2 : NumExp ))
local TI=${TI:-4}
local TF=${TF:-18}
local TI2=${TI2:-$TI}
local TF2=${TF2:-$TF}
local MaxTF=$(( TF2 > TF ? TF2 : TF ))
# Stat - which statistic
# FitOptions - extra options to MultiFit

#LabelTF=${LabelTF:-$((THalf-2))}
local LabelTF=${LabelTF:-$(( MaxTF + 2 ))}

local Corr=${Corr:-corr}

#Stat=${Stat:-0.05} # Doesn't need a default, but used if present

########################################
# Derived

local OutSubDir=2ptp2
local Suffix=$Seed.txt
local Fold=fold.$Suffix
local CorrPrefix=$Corr/2ptp2/${Meson}_p2_${p}_

local FitType=corr_${TI}_${TF}_${TI2}_${TF2}
local MesonPrefix=$Meson/${Meson}_p2_${p}.$FitType.g5P_g5W.model

local CorrFiles="$PlotData/${CorrPrefix}g5P_g5P.$Fold $PlotData/${CorrPrefix}g5P_g5W.$Fold"
local ExtraFiles="$PlotData/${CorrPrefix}gT5P_g5P.$Fold $PlotData/${CorrPrefix}gT5P_g5W.$Fold"

LabelTF="${LabelTF} ${LabelTF}"

local OutFile=$OutSubDir/$MesonPrefix.$Seed
local LogFile=$OutFile.log

########################################
# Perform the fit

local Cmd="MultiFit -e $NumExp --mindp 1"
[ -v Stat ] && Cmd="$Cmd --Hotelling $Stat"
Cmd="$Cmd --overwrite"
#Cmd="$Cmd --debug-signals"
[ -v FitOptions ] && Cmd="$Cmd $FitOptions"
Cmd="$Cmd --summary 2 -i $PlotData/${CorrPrefix} -o $OutSubDir/$Meson/"
Cmd="$Cmd g5P_g5P.fold.$Seed.h5,t=${TI}:${TF},e=$NumExp"
Cmd="$Cmd g5P_g5W.fold.$Seed.h5,t=${TI2}:${TF2},e=$NumExp2"
mkdir -p "$OutSubDir/$Meson"
#echo "$Cmd"
echo "$Cmd"  > $LogFile
if  ! $Cmd &>> $LogFile
then
  LastError=${PIPESTATUS[0]}
  if (( LastError != 3 ))
  then
    echo "MultiFit error $LastError"
    return 1
  fi
  echo "Warning: Not all fit parameters resolved"
fi

########################################
# Plot the fit

# Get the energy difference
GetColumnValues $OutFile.h5 "${Meson//_/-} (n^2=$p) E_0=" '' E0
if ! ColumnValues=$(GetColumn --exact ChiSqPerDof,pValueH --partial E0 $OutFile.h5)
then
  LastError=${PIPESTATUS[0]}
  echo "Error $LastError: $ColumnValues"
  unset ColumnValues
else
  #echo "OK: $ColumnValues"
  ColumnValues=($ColumnValues)
  E0="${ColumnValues[@]:16:8}"
  RefText="${Meson//_/-} (n^2=$p) E_0=${ColumnValues[@]:17:1}, χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})"
fi

# Plot it

Cmd="title='point-point point-wall' files='$CorrFiles' tf='$LabelTF' save='$OutFile'"
Cmd="$Cmd yrange='${ayRange[$Meson,$p]}'"
[ -v RefText ] && Cmd="$Cmd RefText='$RefText' RefVal='${ColumnValues[@]:16:8}'"
[ -v ExtraFiles ] && Cmd="$Cmd extra='$ExtraFiles'"
Cmd="$Cmd plottd.sh $OutSubDir/${MesonPrefix}_td.$Suffix"
DoCmd
}

############################################################

# Scan a two-point fit range

############################################################

function TwoPointScan()
{
  local OptP="$1" # t=...[,e=...]
  local OptW="$2"

  local qHeavy=${qHeavy:-h$Heavy}
  local Meson=${Meson:-${qHeavy}_s}
  local p=${p:-0}

  # SubDir - which subdirectory to place output in
  # SummaryOptions - extra options for the Summary

  local Corr=${Corr:-corr}

  ########################################
  # Error check

  local DoErrorMsg=0
  if [ -z "$OptP" ]
  then
    echo "Point options unspecified. Expect as a minimum t=ti:tf"
    DoErrorMsg=1
  fi
  if [ -z "$OptW" ]
  then
    echo "Wall options unspecified. Expect as a minimum t=ti:tf"
    DoErrorMsg=1
  fi
  if (( DoErrorMsg )); then return 1; fi

  ########################################
  # Derived

  if [[ ! $p =~ _ ]]; then p="p2_$p"; fi # Assume momentum is p^2 (but allow it to be set explicitly)

  local InDir=$PlotData/$Corr/2ptp2
  local OutDir=$Ensemble/TwoPointScan/$Meson
  if [ -v SubDir ]; then OutDir="$OutDir/$SubDir"; fi

  local Suffix=fold.$Seed.h5
  local FileBase=${Meson}_${p}_g5P_g5
  local LogFile=$OutDir/${FileBase}W.log

  ########################################
  # Perform the fit

  echo "$Ensemble: Scan two point ${Meson}_${p}, Point $OptP, Wall $OptW"

  mkdir -p $OutDir # So I can redirect output here

  local Cmd="MultiFit -o \"$OutDir/\" -i \"$InDir/$FileBase\""
  Cmd="$Cmd --debug-signals --mindp 2 --iter 100000 --strict 3"
  [ -v FitOptions ] && Cmd="$Cmd $FitOptions"
  Cmd="$Cmd P.$Suffix,$OptP"
  Cmd="$Cmd W.$Suffix,$OptW"
  DoCmd

  Cmd="FitSummary --all --strict 0 -i \"$OutDir/\" -o \"$OutDir/\"" # --errdig 2 --tablen
  [ -v SummaryOptions ] && Cmd="$Cmd $SummaryOptions"
  Cmd="$Cmd '${Meson}_${p}.*.h5'"
  DoCmd

  Cmd="rm '$OutDir/${Meson}_${p}.corr_'*"
  DoCmd
}

############################################################

# Simultaneous fits of point-point data at all momenta

############################################################

function SimulP()
{
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  sTimes="${aTimes[*]}"
  sTimes=${sTimes// /_}
  sTimes=${sTimes//:/_}
  Cmd="$MultiFit -e 2 -N 24 -o $OutDir/s_l"
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do
    Cmd="$Cmd ${aFitFiles[i]}.h5,t=${aTimes[i]}"
  done
  BaseFile="$OutDir/s_l.corr_$sTimes.g5P.model"
  LogFile="$BaseFile.$Seed.log"
  DoCmd

  # Get the energy difference
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do
    EKeys=$EKeys${EKeys:+,}s_l_p2_${i}-E0
  done
  GetColumnValues $BaseFile.$Seed.h5 "E_0(n^2=0)=" $EKeys

  # Now plot it
  for((i=0; i < ${#aFitFiles[@]}; ++i))
  do
    ThisTF=${aTimes[i]}; ThisTF=${ThisTF##*:}; ThisTF=$((ThisTF+3))
    tf="$tf${tf:+ }$ThisTF"
    title="$title${title:+ }n^2=$i"
    files="$files${files:+ }${aFitFiles[i]}.txt"
  done
  Cmd="tf='$tf' title='$title' files='$files'"
  [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
  Cmd="$Cmd save='$BaseFile.$Seed' yrange=0.3:0.7 plottd.sh ${BaseFile}_td.$Seed.txt"
  DoCmd
  for((i=0; i < ${#aFitFiles[@]}; ++i))
  do
    tf=${aTimes[i]}; tf=${tf##*:}; tf=$((tf+3))
    title="n^2=$i"
    files="${aFitFiles[i]}.txt"
    Cmd="title='$title' files='$files'"
    #Cmd="$Cmd tf='$tf'"
    [ -v RefText ] && Cmd="$Cmd RefText='E_0(n^2=${i})=${ColumnValues[@]:$((17+i*8)):1}, χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})'"
    Cmd="$Cmd save='$OutDir/s_l_p2_${i}.corr_$sTimes.g5P.model.$Seed'"
    Cmd="$Cmd mmin=$i mmax=$i plottd.sh ${BaseFile}_td.$Seed.txt"
    DoCmd
  done
}

############################################################

# For each non-zero momentum
# fit pp & pw data using E(n^2=0) and lattice dispersion relation

# Input:
#  aTimes:  array of point-point start:stop fit times for each momenta (0'th ignored)
#  aTimesW: array of point-wall  start:stop fit times for each momenta (0'th ignored)
#  PriorFitTimes: fit to use for zero momentum

############################################################

function FitEachMomPW()
{
  local OutGroup=$1
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local ZeroMomFit="$OutDir/s_l_p2_0.corr_$PriorFitTimes.g5P_g5W.model.1835672416.h5"
  for((i=1; i < 5; ++i))
  do
    sTimes="${aTimes[i]} ${aTimesW[i]}"
    sTimes=${sTimes// /_}
    sTimes=${sTimes//:/_}
    BaseFile="$OutDir/s_l_p2_$i.${OutGroup}_$PriorFitTimes"
    Cmd="$MultiFit -e 2 -N 24 -o $BaseFile $ZeroMomFit"
    ThisFile="${aFitFiles[i]}.h5"
    Cmd="$Cmd $ThisFile,t=${aTimes[i]}"
    Cmd="$Cmd ${ThisFile//g5P_g5P/g5P_g5W},t=${aTimesW[i]},e=1"
    BaseFile="$BaseFile.corr_$sTimes.g5P_g5W.model"
    LogFile="$BaseFile.$Seed.log"
    DoCmd

    GetColumnValues $BaseFile.$Seed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    ThisFile="${aFitFiles[i]}.txt"
    Cmd="title='point-point point-wall' files='$ThisFile ${ThisFile//g5P_g5P/g5P_g5W}'"
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd save='$BaseFile.$Seed'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$Seed.txt"
    DoCmd
  done
}

############################################################

# For each non-zero momentum
# fit pp & pw data using E(n^2=0) and lattice dispersion relation

# Input:
#  aTimes:  array of point-point start:stop fit times for each momenta (0'th ignored)
#  PriorFitTimes: fit to use for zero momentum

############################################################

function FitEachMomP()
{
  local OutGroup=$1
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local ZeroMomFit="$OutDir/s_l_p2_0.corr_$PriorFitTimes.g5P_g5W.model.1835672416.h5"
  for((i=1; i < 5; ++i))
  do
    sTimes="${aTimes[i]}"
    sTimes=${sTimes//:/_}
    BaseFile="$OutDir/s_l_p2_$i.${OutGroup}_$PriorFitTimes"
    Cmd="$MultiFit -e 2 -N 24 -o $BaseFile $ZeroMomFit"
    Cmd="$Cmd ${aFitFiles[i]}.h5,t=${aTimes[i]}"
    BaseFile="$BaseFile.corr_$sTimes.g5P_g5W.model"
    LogFile="$BaseFile.$Seed.log"
    DoCmd

    GetColumnValues $BaseFile.$Seed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    Cmd="title='point-point n^2=$i' files='${aFitFiles[i]}.txt'"
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd save='$BaseFile.$Seed'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$Seed.txt"
    DoCmd
  done
}
