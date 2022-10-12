#!/usr/bin/env bash

# Prototype used on M1

ModelDir=fit/2ptp2
Models=$ModelDir/h447_s/h447_s_p2_0.corr_9_25_9_25.g5P_g5W.model.1835672416.h5
Models="$Models $ModelDir/s_l/s_l_p2_0.corr_6_23_6_23.g5P_g5W.model.1835672416.h5"

CorrPrefix=corr/3sm_sp2/quark_l_h447_gT_dt_
CorrSuffix=_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5,snk=s_l_p2_0,src=h447_s_p2_0

dt=(20 24 28 32)
MinDP=2
tStart=2 # How far in from each edge
tDelta=2 # How close must each other range be

############################################################

# Parameters end

############################################################

dtLow=${dt[0]}
dtHigh=${dt[-1]}

echo "From $dtLow to $dtHigh"

# This is the range agreed with Tobi
#range=(R-1:-2:2:5:5 R-2:-2:6:5:5 R-3:-2:10:5:5 R-4:-2:14:5:5)
function FitRangeTobi
{
    unset List
    for (( ti=tStart; ti <= dtLow - tStart - MinDP + 1; ++ti ))
    do
	for (( tf=ti + MinDP - 1; tf <= dtLow - tStart - MinDP + 1; ++tf ))
	do
	    List="${List:+${List} }$ti:$tf"
	done
    done
}

# This is all fit ranges
#range=(2:4:15:15 2:4:19:19 2:4:23:23 2:4:27:27)

FitRangeTobi
dtidtf=$((2 * tDelta + 1)):$((2 * tDelta + 1))
for range in $List
do
  tNums=(${range//:/ })
  ti=${tNums[0]}
  tf=${tNums[1]}
  Cmd="MultiFit -i ../ -e 2 --mindp $MinDP $Models"
  Cmd="$Cmd ${CorrPrefix}${dtLow}${CorrSuffix},t=${range}"
  for (( i=1; i < ${#dt[@]}; ++i ))
  do
    RelRange=-${tDelta}:$((dt[i] - dtLow - tDelta)):$dtidtf
    Cmd="$Cmd ${CorrPrefix}${dt[i]}${CorrSuffix},t=R-${i}:$RelRange"
  done
  echo $Cmd
done
