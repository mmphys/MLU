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
RemoteBase=/tessfs1/work/dp008/dp008/shared/data/semilep/${BaseDir}_64
LocalBase=~/data/Study3/$BaseDir
mkdir -p $LocalBase
cd $LocalBase

function GetFiles()
{
  local Dir=$1
  local What=$2
  mkdir -p $Dir
  rsync -uvlmrt tess:"${RemoteBase}/$Dir/$What" $Dir/
}

function GetDebug()
{
  rsync -uvlmrt tess:"${RemoteBase}/corr/2pt/*_p_[01]_0_0_g5P_*.h5" corr/2pt/
  rsync -uvlmrt tess:"${RemoteBase}/corr/3pt_s/*_ps_0_0_0_gT_dt_*_p_0_0_0_g5P_g5P.*.h5 ${RemoteBase}/corr/3pt_s/*_ps_*_0_0_g[TX]_dt_*_p_1_0_0_g5P_g5P.*.h5 ${RemoteBase}/corr/3pt_s/*_ps_-1_0_0_g[TX]_dt_*_p_0_0_0_g5P_g5P.*.h5" corr/3pt_s/
  #rsync -uvlmrt tess:"${RemoteBase}/fit/2pt/*_p_[01]_0_0.corr_* ${RemoteBase}/fit/2pt/*_p_[01]_0_0.corr.*params*" fit/2pt/
}

echo "Getting correlators"
mkdir -p corr
rsync -uvlmrt --include='*.txt' --include='h447_s_p2_0_g5P_g5[PW].*.h5' \
  --include='h447_l_p2_0_g5P_g5[PW].*.h5' \
  --include='l_s_p2_0_g5P_g5[PW].*.h5' \
  --include='l_l_p2_0_g5P_g5[PW].*.h5' \
  --include='quark_h447_h447_gT_dt_20_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5' \
  --include='quark_h447_l_gT_dt_20_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5' \
  --include='quark_l_h447_gT_dt_20_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5' \
  --include='quark_l_l_gT_dt_20_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5' \
  --include='quark_s_s_gT_dt_20_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5' \
  --include='*/' --exclude='*' tess:"${RemoteBase}/corr/" corr/

echo "Getting fits"
mkdir -p fit
rsync -uvlmrt --include='*.params*.txt' --include='*/' \
      --include='h447_l_p2_0.corr_10_20.g5P_g5W.model.*.h5' \
      --include='l_s_p2_0.corr_7_18.g5P_g5W.model.*.h5' \
      --include='l_l_p2_0.corr_6_15.g5P_g5W.model.*.h5' \
      --include='l_l_p2_1.corr_6_15.g5P_g5W.model.*.h5' \
      --exclude='*' tess:"${RemoteBase}/fit/2ptp2/" fit/2ptp2/

mkdir -p ratio/3pt_s/raw
# GetDebug

