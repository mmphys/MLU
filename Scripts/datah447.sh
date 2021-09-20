#!/usr/bin/env bash

#set -x # Log all steps - debug
set -e # Fail on all errors
shopt -s nullglob # Return nothing if no match

if (( $# > 1 )); then
  echo "Copy h447 data from Tesseract - no parameters expected"
  exit
fi

# Set up Remote and Local directories
BaseDir=M1/analyse/h447
RemoteBase=/tessfs1/work/dp008/dp008/shared/data/semilep/${BaseDir}
LocalBase=~/data/Study3/$BaseDir
mkdir -p $LocalBase
cd $LocalBase

#function GetFiles()
#{
#  local Dir=$1
#  local What=$2
#  mkdir -p $Dir
#  rsync -uvlmrt tess:"${RemoteBase}/$Dir/$What" $Dir/
#}

function GetDebugWard()
{
echo "Getting Debug info:"
echo "    - Ward identity bootstraps"
ThisDir=bootstrap
mkdir -p $ThisDir
rsync -uvlmrt --include='PJ5q_*.txt' --include='PA0_*.txt' \
  --include='Ward/' --exclude='*' tess:"${RemoteBase}/$ThisDir/" $ThisDir/
}

function GetDebugMeson()
{
echo "Getting Debug info:"
echo "    - Mesons 3pt"
local ThisDir=~/data/Study3/M1/M1.cpu/meson/3pt_s
mkdir -p $ThisDir
rsync -uvlmrt --include='PP_quark_rev_l_h447_g5_dt_16_ps_-*_p_0_0_0_t_48.*.h5' --exclude='*' \
  tess:"/tessfs1/work/dp008/dp008/shared/data/semilep/M1/M1.cpu/900/3pt_s/ZFPP_quark_rev/" $ThisDir/
}

function GetCorrelators()
{
echo "Getting correlators"
ThisDir=corr
mkdir -p $ThisDir
rsync -uvlmrt --include='*.txt' \
  --include='PJ5q_*.h5' \
  --include='h447_h447_p_0_0_0_g5[PW]_g5[PW].*.h5' \
  --include='l_l_p_0_0_0_g5[PW]_g5[PW].*.h5' \
  --include='s_s_p_0_0_0_g5[PW]_g5[PW].*.h5' \
  --include='h447_h447_p2_0_g5[PW]_g5[PW].*.h5' \
  --include='[ls]_[ls]_p2_0_g5[PW]_g5[PW].*.h5' \
  --include='h447_[ls]_p2_[01]_g5[PW]_g5[PW].*.h5' \
  --include='l_[ls]_p2_[01]_g5[PW]_g5[PW].*.h5' \
  --include='quark_h447_h447_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_h447_l_gT_dt_20_p2_[01]_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_l_h447_gT_dt_20_p2_0_ps2_[01]_g5*_g5*.fold.1835672416.h5' \
  --include='quark_l_l_gT_dt_20_p2_[01]_ps2_[01]_g5*_g5*.fold.1835672416.h5' \
  --include='quark_s_s_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='*/' --exclude='*' tess:"${RemoteBase}/$ThisDir/" $ThisDir/
}

function GetFits()
{
echo "Getting fits"
ThisDir=fit
mkdir -p $ThisDir
rsync -uvlmrt --exclude='bad*/' --include='*/' --include='*.params*.txt' --include='Fit*.txt' \
      --include='ZV_*.bootstrap.*.txt' \
      --include='h447_l_p2_0.corr_10_16_11_20.g5P_g5W.model.1835672416.h5' \
      --include='h447_s_p2_0.corr_10_19_12_19.g5P_g5W.model.1835672416.h5' \
      --include='l_s_p2_*.corr_6_11_4_9.g5P_g5W.model.1835672416.h5' \
      --include='l_l_p2_*.corr_6_11_3_8.g5P_g5W.model.1835672416.h5' \
      --include='ZV_l_dt_28.corr_10_18.g5P_g5P.model.1835672416.h5' \
      --include='ZV_h447_dt_24.corr_10_14.g5P_g5P.model.1835672416.h5' \
      --exclude='*' tess:"${RemoteBase}/$ThisDir/" $ThisDir/
}

function GetRatios()
{
echo "Getting ratios"
ThisDir=ratio
mkdir -p $ThisDir
rsync -uvlmrt --include='*.txt' --include='*/' --exclude='*' tess:"${RemoteBase}/$ThisDir/" $ThisDir/
}

#GetDebugWard
#GetDebugMeson
GetCorrelators
#GetFits # These can be a little slow when there's 50,000 of them
#GetRatios
