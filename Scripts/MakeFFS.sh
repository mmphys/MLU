#!/bin/env bash

#set -x
set -e

Analyse=analyse

function EchoCmd()
{
  echo $Cmd
  echo $Cmd &>> log/MakeFFS.log
}

function DoCmd()
{
  EchoCmd
  eval $Cmd &>> log/MakeFFS.log
}

if ! [ -d "$Analyse" ]
then
  echo "$Analyse directory missing"
  exit 1
fi

Ensemble=${PWD##*/}
unset N
case $Ensemble in
  C1 | C2) N=24;;
  M1 | M2 | M3) N=32;;
  F1M) N=48;;
esac

if ! [ -v N ]; then
  echo "Ensemble $Ensemble unrecognised"
  exit 1
fi

mkdir -p "$Analyse/ffs/log"
cd "$Analyse/ffs"
Cmd="Making form factors for $Ensemble, N=$N"; EchoCmd
for Dir in lp2 sp2
do
  Cmd="CRatio --i2 ../corr/2ptp2/ --i3 '../RatioFit/$Dir/' -o '$Dir/' --efit '../fit/Fit_$Dir.txt' --type f,$N '*_gT_*.h5'"
  DoCmd
  Cmd="FitSummary -o '$Dir/Summary/' '$Dir/*.h5'"
  DoCmd
done
