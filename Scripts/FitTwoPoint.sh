#!/usr/bin/env bash

# Simple two-point plots

set -e

############################################################

# Inputs

############################################################

#Ensemble=${Ensemble:-F1M}
. PlotCommon.sh
#set -x

qHeavy=${qHeavy:-h$Heavy}
Meson=${Meson:-${qHeavy}_s}
p=${p:-0}
NumExp=${NumExp:-3}
NumExp2=${NumExp2:-${NumExp}}
MaxExp=$(( NumExp2 > NumExp ? NumExp2 : NumExp ))
TI=${TI:-4}
TF=${TF:-18}
TI2=${TI2:-$TI}
TF2=${TF2:-$TF}
MaxTF=$(( TF2 > TF ? TF2 : TF ))
# Stat - which statistic
# FitOptions - extra options to MultiFit

#LabelTF=${LabelTF:-$((THalf-2))}
LabelTF=${LabelTF:-$(( MaxTF + 2 ))}

Seed=${Seed:-1835672416}
Corr=${Corr:-corr}
MELFit=${MELFit:-MELFit}

#Stat=${Stat:-0.05} # Doesn't need a default, but used if present

############################################################

# Derived

############################################################

OutSubDir=$Ensemble/$MELFit/2ptp2

Suffix=$Seed.txt
Fold=fold.$Suffix
CorrPrefix=$Corr/2ptp2/${Meson}_p2_${p}_

FitType=corr_${TI}_${TF}_${TI2}_${TF2}
MesonPrefix=$Meson/${Meson}_p2_${p}.$FitType.g5P_g5W.model

CorrFiles="$PlotData/${CorrPrefix}g5P_g5P.$Fold $PlotData/${CorrPrefix}g5P_g5W.$Fold"
ExtraFiles="$PlotData/${CorrPrefix}gT5P_g5P.$Fold $PlotData/${CorrPrefix}gT5P_g5W.$Fold"

LabelTF="${LabelTF} ${LabelTF}"

OutFile=$OutSubDir/$MesonPrefix.$Seed

############################################################

# Perform the fit

############################################################

Cmd="MultiFit -e $NumExp --mindp 1"
[ -v Stat ] && Cmd="$Cmd --Hotelling $Stat"
Cmd="$Cmd --overwrite"
#Cmd="$Cmd --debug-signals"
[ -v FitOptions ] && Cmd="$Cmd $FitOptions"
Cmd="$Cmd --summary 2 -i $PlotData/${CorrPrefix} -o $OutSubDir/$Meson/"
Cmd="$Cmd g5P_g5P.fold.$Seed.h5,t=${TI}:${TF},e=$NumExp"
Cmd="$Cmd g5P_g5W.fold.$Seed.h5,t=${TI2}:${TF2},e=$NumExp2"
mkdir -p "$OutSubDir/$Meson"
#echo "$Cmd"
echo "$Cmd"  > $OutFile.log
if  ! $Cmd &>> $OutFile.log
then
  LastError=${PIPESTATUS[0]}
  if (( LastError != 3 ))
  then
    echo "MultiFit error $LastError"
    exit 1
  fi
  echo "Warning: Not all fit parameters resolved"
fi

############################################################

# Plot the fit

############################################################

# Get the energy difference
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
[ -v yrange ] && Cmd="$Cmd yrange='$yrange'"
[ -v E0 ] && Cmd="$Cmd RefVal='$E0'"
[ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
[ -v ExtraFiles ] && Cmd="$Cmd extra='$ExtraFiles'"
Cmd="$Cmd plottd.sh $OutSubDir/${MesonPrefix}_td.$Suffix"
#echo "$Cmd"
echo "$Cmd" >> $OutFile.log
eval  $Cmd
