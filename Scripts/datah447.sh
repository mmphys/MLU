#!/usr/bin/env bash

set -x # Log all steps - debug
set -e # Fail on all errors
shopt -s nullglob # Return nothing if no match

if (( $# > 1 )); then
  echo "Copy h447 data from Tesseract - no parameters expected"
  exit
fi

BaseDir=data/Study3/M1/analyse/h447
mkdir -p ~/$BaseDir
cd ~/$BaseDir

mkdir -p corr/2pt corr/3pt_s fit/2pt ratio/3pt_s/raw

rsync -uvlmrt tess:"${BaseDir}/corr/2pt/*_p_[01]_0_0_g5P_*" corr/2pt/
rsync -uvlmrt tess:"${BaseDir}/corr/3pt_s/*_ps_0_0_0_gT_dt_*_p_0_0_0_g5P_g5P.* ${BaseDir}/corr/3pt_s/*_ps_*_0_0_g[TX]_dt_*_p_1_0_0_g5P_g5P.* ${BaseDir}/corr/3pt_s/*_ps_-1_0_0_g[TX]_dt_*_p_0_0_0_g5P_g5P.*" corr/3pt_s/
rsync -uvlmrt tess:"${BaseDir}/fit/2pt/*_p_[01]_0_0.corr_* ${BaseDir}/fit/2pt/*_p_[01]_0_0.corr.*params*" fit/2pt/
