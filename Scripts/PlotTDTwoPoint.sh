#!/usr/bin/env bash

# Simple two-point plots
. PlotCommon.sh
#set -x
set -e

############################################################

# Inputs

############################################################

Meson=${Meson:-s_l}
p=${p:-0}
Ensemble=${Ensemble:-F1M}
InBase=${InBase:-/Volumes/QCD/tursa/semilep/data}
NumExp=${NumExp:-2}
TI=${TI:-7}
TF=${TF:-18}
LabelTF=${LabelTF:-30}

Seed=${Seed:-1835672416}
Analyse=${Analyse:-analyse}
Corr=${Corr:-corr}
Fit=${Fit:-fit}
Plot=${Plot:-Plot}

############################################################

# Derived

############################################################

if [ "$NumExp" != 2 ]; then Fit=$Fit$NumExp; fi

FitType=corr_${TI}_${TF}_$((TI-1))_${TF}

DataDir=$InBase/$Ensemble/$Analyse

Suffix=$Seed.txt
CorrPrefix=$Corr/2ptp2/${Meson}_p2_${p}_
MesonPrefix=$Fit/2ptp2/$Meson/${Meson}_p2_${p}.$FitType.g5P_g5W.model

CorrFiles="$DataDir/${CorrPrefix}g5P_g5P.fold.$Suffix $DataDir/${CorrPrefix}g5P_g5W.fold.$Suffix"
ExtraFiles="$DataDir/${CorrPrefix}gT5P_g5P.fold.$Suffix" #$DataDir/${CorrPrefix}g5P_g5W.fold.$Suffix"

LabelTF="${LabelTF} ${LabelTF}"

OutSubDir=$Plot/$Fit/$Ensemble
OutFile=$OutSubDir/${Meson}_p2_${p}.E${NumExp}_${FitType}

############################################################

# Plot the fit

############################################################

# Get the energy difference
ColumnValues=$(GetColumn --exact ChiSqPerDof --partial E0 $DataDir/$MesonPrefix.$Seed.h5)
if [ "$?" != 0 ]; then
  echo "Error: $ColumnValues"
  unset ColumnValues
else
  echo "OK: $ColumnValues"
  ColumnValues=($ColumnValues)
  ChiSqPerDof="χ²/dof=${ColumnValues[@]:3:1}"
  E0="${ColumnValues[@]:7:7}"
fi

# Plot it

mkdir -p $OutSubDir
Cmd="title='point-point point-wall' files='$CorrFiles' tf='$LabelTF' save='$OutFile'"
[ -v yrange ] && Cmd="$Cmd yrange='$yrange'"
[ -v E0 ] && Cmd="$Cmd RefVal='$E0'"
[ -v ChiSqPerDof ] && Cmd="$Cmd RefText='$ChiSqPerDof'"
[ -v ExtraFiles ] && Cmd="$Cmd extra='$ExtraFiles'"
Cmd="$Cmd plottd.sh $DataDir/${MesonPrefix}_td.$Suffix"
#echo "$Cmd"
eval  $Cmd
