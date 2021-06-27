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

function GetFiles()
{
  local Dir=$1
  local What=$2
  mkdir -p $Dir
  rsync -uvlmrt tess:"${RemoteBase}/$Dir/$What" $Dir/
}

function GetDebug()
{
echo "Getting Debug info:"
echo "    - Ward identity bootstraps"
rsync -uvlmrt --include='PJ5q_*.txt' --include='PA0_*.txt' \
  --include='Ward/' --exclude='*' tess:"${RemoteBase}/bootstrap/" bootstrap/
}

echo "Getting correlators"
mkdir -p corr
rsync -uvlmrt --include='*.txt' \
  --include='h447_[ls]_p2_[01]_g5[PW]_g5[PW].*.h5' \
  --include='l_[ls]_p2_0_g5[PW]_g5[PW].*.h5' \
  --include='quark_h447_h447_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_h447_l_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_l_h447_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_l_l_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='quark_s_s_gT_dt_20_p2_0_ps2_0_g5*_g5*.fold.1835672416.h5' \
  --include='*/' --exclude='*' tess:"${RemoteBase}/corr/" corr/

echo "Getting fits"
mkdir -p fit
rsync -uvlmrt --include='*.params*.txt' --include='*/' \
      --include='h447_l_p2_0.corr_10_20.g5P_g5W.model.*.h5' \
      --include='l_s_p2_0.corr_7_18.g5P_g5W.model.*.h5' \
      --include='l_l_p2_0.corr_6_15.g5P_g5W.model.*.h5' \
      --include='l_l_p2_1.corr_6_15.g5P_g5W.model.*.h5' \
      --exclude='*' tess:"${RemoteBase}/fit/2ptp2/" fit/2ptp2/

mkdir -p ratio/3pt_s/raw bootstrap
GetDebug

