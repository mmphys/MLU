#!/usr/bin/env bash

############################################################

# Make ratios in PlotData using fits in NoSync directories

# Input:
#   Ensemble: The ensemble to make ratios for
#   series:   (default: renorm) which two-pt fits enter the extraction of Z_V
#   ZVEns:    (optional) the ensemble the ZVfits come from
#   Suffix:   (optional) Suffix for: 1) ZV$Suffix.txt; 2) ratio$Suffix

############################################################

#set -x
set -e

# Validate parameters

series=${series:-renorm}
BaseDir="$PWD"
EnsembleDir="$BaseDir/$Ensemble"
ZVEnsRoot="$BaseDir/${ZVEns:-$Ensemble}"
ZV="$ZVEnsRoot/Renorm/ZV$Suffix.txt"
EFit="$ZVEnsRoot/MELFit/Fit_sp2_$series.txt"

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

SubDir=3sm_sp2
Corr2="$PlotData/corr/2ptp2"
Corr3="$PlotData/corr/$SubDir"
Ratio="$PlotData/ratio$Suffix/$SubDir"
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
