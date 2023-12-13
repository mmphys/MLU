#!/usr/bin/env bash

############################################################

# Make ratios in PlotData using fits in NoSync directories

# Input:
#   Ensemble: The ensemble to make ratios for
#   ZV:       (optional) the ensemble the ZVfits come from

############################################################

#set -x
set -e

# Validate parameters

series=${series:-renorm}
BaseDir="$PWD"
EnsembleDir="$BaseDir/$Ensemble"
ZV="$BaseDir/${ZV:-$Ensemble}/Renorm/ZV.txt"

if ! [ -v Ensemble ]; then
  echo "Ensemble not set"
  exit 2
fi
if ! [ -d "$EnsembleDir" ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit 2
fi

if ! [ -e "$ZV" ]; then
  echo "ZV $ZV doesn't exist"
  exit 2
fi

############################################################

# Run

############################################################

. PlotCommon.sh

EFit="$EnsembleDir/MELFit/Fit_sp2_$series.txt"
SubDir=3sm_sp2
Corr2="$PlotData/corr/2ptp2"
Corr3="$PlotData/corr/$SubDir"
Ratio="$PlotData/ratio/$SubDir"
LogFile="$Ratio/MakeRatio.log"

if ! [ -d "$EFit" ] && ! [ -h "$EFit" ]; then
  echo "$EFit doesn't exist"
  exit 2
fi
if ! [ -d "$PlotData" ]; then
  echo "$PlotData doesn't exist"
  exit 2
fi
if [ -d "$Ratio" ]; then
  echo "$Ratio already exists"
  exit 2
fi

mkdir -p "$Ratio"

Cmd="CRatio --i2 '$Corr2/' --i3 '$Corr3/' -o '$Ratio/' --efit '$EFit' --type R,$ZV"
Cmd+=" '*_l_h*_g*_dt_*_p2_0_ps2_*_g5P_g5P.*.h5'"

{
  echo ZV=$ZV
  cat "$ZV"
  echo EFit=$EFit
  cat "$EFit"
  echo "$Cmd"
  eval "$Cmd"
} > $LogFile
