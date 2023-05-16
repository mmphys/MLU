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
  if [[ $LogFile != ${LogFile%/*} && -n ${LogFile%/*} ]]; then mkdir -p ${LogFile%/*}; fi
  if [ "$1" = EraseLog ] && [ -f $LogFile ]; then rm $LogFile; fi
  #echo $Cmd
  echo $Cmd &>> $LogFile
  local RetVal=0
  if ! eval $Cmd &>> $LogFile; then
    RetVal=${PIPESTATUS[0]}
    # Ignore MultiFit 'Not all params resolved" warning
    if [[ ${Cmd:0:8} = MultiFit && $RetVal = 3 ]]; then
      RetVal=0
      echo "Warning: not all parameters resolved: $Cmd"
    else
      echo "Error $RetVal" &>> $LogFile
      echo "Error $RetVal executing $Cmd"
    fi
  fi
  return $RetVal
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
local Thin1=${Thin1+,t$Thin1}
local Thin2=${Thin2+,t$Thin2}

#LabelTF=${LabelTF:-$((THalf-2))}
local LabelTF=${LabelTF:-$(( MaxTF + 2 ))}

local Corr=${Corr:-corr}

#Stat=${Stat:-0.05} # Doesn't need a default, but used if present

########################################
# Derived

local OutSubDir=$Ensemble/MELFit/2ptp2/$Meson
local Fold=fold.$DataSeed
local CorrPrefix=$Corr/2ptp2/${Meson}_p2_${p}_

local FitType=corr_${TI}_${TF}_${TI2}_${TF2}
local MesonPrefix=${Meson}_p2_${p}.$FitType.g5P_g5W.model

local ExtraFiles="$PlotData/${CorrPrefix}gT5P_g5P.$Fold.txt $PlotData/${CorrPrefix}gT5P_g5W.$Fold.txt"

LabelTF="${LabelTF} ${LabelTF}"

local OutFile=$OutSubDir/$MesonPrefix.$MLUSeed
local LogFile=$OutFile.log
mkdir -p "$OutSubDir"

########################################
# Perform the fit

local Cmd="MultiFit -e $NumExp --mindp 1"
[ -v Stat ] && Cmd="$Cmd --Hotelling $Stat"
Cmd="$Cmd --overwrite"
Cmd="$Cmd --debug-signals"
Cmd="$Cmd --strict 1"
[ -v FitOptions ] && Cmd="$Cmd $FitOptions"
Cmd="$Cmd --summary 2 -i $PlotData/${CorrPrefix} -o $OutSubDir/"
Cmd="$Cmd g5P_g5P.$Fold.h5,t=${TI}:${TF}${Thin1},e=$NumExp"
Cmd="$Cmd g5P_g5W.$Fold.h5,t=${TI2}:${TF2}${Thin2},e=$NumExp2"
#echo "$Cmd"
DoCmd EraseLog || return 0

########################################
# Plot the fit

# Get E0 (reference value for plot)
GetColumnValues $OutFile.h5 "${Meson//_/-} (n^2=$p) E_0=" '' E0
Cmd="title='point-point point-wall' tf='$LabelTF' yrange='${ayRange[$Meson,$p]}'"
[ -v RefText ] && Cmd="$Cmd RefText='$RefText' RefVal='${ColumnValues[@]:16:8}'"
[ -v ExtraFiles ] && Cmd="$Cmd extra='$ExtraFiles'"
Cmd="$Cmd plottd.sh $OutSubDir/${MesonPrefix}_td.$MLUSeed.txt"
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

  local Suffix=fold.$MLUSeed.h5
  local FileBase=${Meson}_${p}_g5P_g5
  local LogFile=$OutDir/${FileBase}W.log

  ########################################
  # Perform the fit

  echo "$Ensemble: Scan two point ${Meson}_${p}, Point $OptP, Wall $OptW"

  mkdir -p $OutDir # So I can redirect output here

  local Cmd="MultiFit -o \"$OutDir/\" -i \"$InDir/$FileBase\""
  Cmd="$Cmd --debug-signals --mindp 2 --iter 100000 --strict 3 --covsrc binned"
  [ -v FitOptions ] && Cmd="$Cmd $FitOptions"
  Cmd="$Cmd P.$DataSeed.txt,$OptP"
  Cmd="$Cmd W.$DataSeed.txt,$OptW"
  DoCmd || return 0

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
  local i EKeys title
  sTimes="${aTimes[*]}"
  sTimes=${sTimes// /_}
  sTimes=${sTimes//:/_}
  Cmd="$MultiFit -e 2 -N $L -o $OutDir/s_l"
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do
    Cmd="$Cmd ${aFitFiles[i]}.h5,t=${aTimes[i]}${aThin[i]:+t${aThin[i]}}"
  done
  BaseFile="$OutDir/s_l.corr_$sTimes.g5P.model"
  LogFile="$BaseFile.$MLUSeed.log"
  DoCmd EraseLog || return 0

  # Get the energy difference
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do
    EKeys=$EKeys${EKeys:+,}s_l_p2_${i}-E0
  done
  GetColumnValues $BaseFile.$MLUSeed.h5 "E_0(n^2=0)=" $EKeys

  # Now plot it
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do title="$title${title:+ }n^2=$i"; done
  Cmd="title='$title'"
  [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
  Cmd="$Cmd plottd.sh ${BaseFile}_td.$MLUSeed.txt"
  DoCmd
  for((i=0; i < ${#aFitFiles[@]}; ++i))
  do
    Cmd="title='n^2=$i'"
    [ -v RefText ] && Cmd="$Cmd RefText='E_0(n^2=${i})=${ColumnValues[@]:$((17+i*8)):1}, χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd save='$OutDir/s_l_p2_${i}.corr_$sTimes.g5P.model.$MLUSeed'"
    Cmd="$Cmd mmin=$i mmax=$i plottd.sh ${BaseFile}_td.$MLUSeed.txt"
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
    LogFile="$BaseFile.$MLUSeed.log"
    DoCmd || return 0

    GetColumnValues $BaseFile.$MLUSeed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    Cmd="title='point-point point-wall' "
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$MLUSeed.txt"
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
    LogFile="$BaseFile.$MLUSeed.log"
    DoCmd || return 0

    GetColumnValues $BaseFile.$MLUSeed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    Cmd="title='point-point n^2=$i'"
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$MLUSeed.txt"
    DoCmd
  done
}
