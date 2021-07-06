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
mkdir -p bootstrap
rsync -uvlmrt --include='PJ5q_*.txt' --include='PA0_*.txt' \
  --include='Ward/' --exclude='*' tess:"${RemoteBase}/bootstrap/" bootstrap/
}

function GetDebugMeson()
{
echo "Getting Debug info:"
echo "    - Mesons 3pt"
local ThisDir=../../
mkdir -p bootstrap/3pt
rsync -uvlmrt --include='PP_quark_rev_l_h447_g5_dt_16_ps_-*_p_0_0_0_t_48.*.h5' --exclude='*' \
  tess:"/tessfs1/work/dp008/dp008/shared/data/semilep/M1/M1.cpu/900/3pt_s/ZFPP_quark_rev/" \
  ~/data/Study3/M1/M1.cpu/meson/3pt_s/
}

function GetCorrelators()
{
echo "Getting correlators"
mkdir -p corr
rsync -uvlmrt --include='*.txt' \
  --include='h447_[ls]_p2_[01]_g5[PW]_g5[PW].*.h5' \
  --include='l_[ls]_p2_[01]_g5[PW]_g5[PW].*.h5' \
  --include='quark_h447_h447_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_h447_l_gT_dt_20_p2_[01]_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_l_h447_gT_dt_20_p2_0_ps2_[01]_g5*_g5*.fold.1835672416.h5' \
  --include='quark_l_l_gT_dt_20_p2_[01]_ps2_[01]_g5*_g5*.fold.1835672416.h5' \
  --include='quark_s_s_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='*/' --exclude='*' tess:"${RemoteBase}/corr/" corr/
}

function GetFits()
{
echo "Getting fits"
mkdir -p fit
rsync -uvlmrt --include='*.params*.txt' --include='Fit*.txt' --include='*/' \
      --include='ZV_*_dt_20.*.model.*.h5' \
      --include='h447_l_p2_0.corr_9_20.g5P_g5W.model.*.h5' \
      --include='l_s_p2_0.corr_7_16.g5P_g5W.model.*.h5' \
      --include='l_l_p2_0.corr_5_14.g5P_g5W.model.*.h5' \
      --include='l_l_p2_1.corr_5_15.g5P_g5W.model.*.h5' \
      --exclude='*' tess:"${RemoteBase}/fit/" fit/
}

function GetRatios()
{
echo "Getting ratios"
mkdir -p ratio
rsync -uvlmrt --include='*.txt' --include='*/' --exclude='*' tess:"${RemoteBase}/ratio/" ratio/
}

#GetDebugWard
GetDebugMeson
#GetCorrelators
#GetFits
#GetRatios
