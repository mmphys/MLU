#!/usr/bin/env bash

set -e

############################################################
#
# Generate a list of fit ranges
# This is the range agreed with Tobi
#
# Entry
#   tStart How far in from each end
#   dtLow   Smallest DeltaT
#   dtHigh  Largest DeltaT
#   MinDP   Minimum number of data points in each correlator
#
# Exit
#   List    List of fit ranges for dtLow
#
############################################################

#range=(R-1:-2:2:5:5 R-2:-2:6:5:5 R-3:-2:10:5:5 R-4:-2:14:5:5)
function FitRangeTobi
{
  local ti
  local tf
  unset List
  for (( ti=tStart; ti <= dtLow - tStart - MinDP + 1; ++ti ))
  do
	for (( tf=ti + MinDP - 1; tf <= dtLow - tStart - MinDP + 1; ++tf ))
	do
	    List="${List:+${List} }$ti:$tf"
	done
  done
}

############################################################
#
# Output fit commands for one decay
#
############################################################

function OutFitCommands
{
  local SpecDir=$1; shift
  local NumExp=$1; shift
  local CorrBase=$1; shift
  local Decay=$1; shift
  local Gamma=$1; shift
  local dtString=$1; shift
  local SnkMeson=$1; shift
  local SnkMom=$1; shift
  local SnkFit=$1; shift
  local SrcMeson=$1; shift
  local SrcMom=$1; shift
  local SrcFit=$1; shift

  local dt=($dtString)
  dtString=${dtString// /_}

  local SnkMesonMom=${SnkMeson}_p2_${SnkMom}
  local SrcMesonMom=${SrcMeson}_p2_${SrcMom}

  local BDG=${CorrBase}_${Decay}_${Gamma}
  local Mom3pt=p2_${SrcMom}_ps2_${SnkMom}

  local CorrPrefix=corr/$SpecDir/${BDG}_dt_
  local CorrSuffix=g5P_g5P.fold.$Suffix
        CorrSuffix=_${Mom3pt}_${CorrSuffix},snk=${SnkMesonMom},src=${SrcMesonMom}
  local ModelSuffix=g5P_g5W.model.$Suffix

  local Models=$ModelDir/$SnkMeson/${SnkMesonMom}.corr_${SnkFit}_${SnkFit}.$ModelSuffix
  Models="$Models $ModelDir/$SrcMeson/${SrcMesonMom}.corr_${SrcFit}_${SrcFit}.$ModelSuffix"

  local dtLow=${dt[0]}
  local dtHigh=${dt[-1]}
  local dtidtf=$((2 * tDelta + 1)):$((2 * tDelta + 1))
  local tNums
  local ti
  local tf
  local range
  local Cmd
  local RelRange

  local OutDir=$SpecDir/${BDG}_${Mom3pt}/dt_$dtString/E${NumExp}_${SnkFit}_${SrcFit}
  local OutSub

  FitRangeTobi
  for range in $List
  do
    tNums=(${range//:/ })
    ti=${tNums[0]}
    tf=${tNums[1]}
    local Sub="${ti}_${tf}"

    Cmd="MultiFit -e $NumExp --mindp $MinDP -i ../ -o $OutDir/$Sub/MEL $Models"
    Cmd="$Cmd ${CorrPrefix}${dtLow}${CorrSuffix},t=${range}"
    for (( i=1; i < ${#dt[@]}; ++i ))
    do
      RelRange=-${tDelta}:$((dt[i] - dtLow - tDelta)):$dtidtf
      Cmd="$Cmd ${CorrPrefix}${dt[i]}${CorrSuffix},t=R-${i}:$RelRange"
    done
    echo $Cmd >> $Analyse1
    OutSub="$OutSub${OutSub:+,}${Sub}"
  done
  echo "FitSummary -i $OutDir/ -o $OutDir/ {$OutSub}'*.model.$Suffix'" >> $Analyse2
}

############################################################

# Start here

############################################################

Heavy=h385
ModelDir=fit/2ptp2
Suffix=1835672416.h5

MinDP=2 # Minimum number of data points in each correlator
tStart=2 # How far in from each edge
tDelta=2 # How close must each other range be

Analyse=analyse.MELFit
Analyse1=$Analyse.sh
Analyse2=$Analyse.2.sh
if [ -f $Analyse1 ]; then rm $Analyse1; fi
if [ -f $Analyse2 ]; then rm $Analyse2; fi

#echo "From $dtLow to $dtHigh"

for deltaT in '16 20 24 28 32' #'16 20 24 28' '20 24 28 32'
do
  for SnkFit in 6_19 9_22
  do
    for SrcFit in 9_22 12_25
    do
      OutFitCommands 3sm_sp2 2 quark l_${Heavy} gT "$deltaT" s_l 0 $SnkFit ${Heavy}_s 0 $SrcFit
    done
  done
done
