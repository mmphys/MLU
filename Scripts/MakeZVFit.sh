#!/bin/env bash

# Make an analysis job which will perform fits of the ZV
# Optional environment variables:
JobRoot=${JobRoot:-analyse}
JobIn=${In:-ZV}
JobOut=${Out:-ZVFit}

# How to reach correlators during run
InDir=../$JobIn
# How to reach directories now (while this script runs)
InDirNow=$JobRoot/$JobIn
OutDirNow=$JobRoot/$JobOut

# Don't expect these to be overriden
JobBaseName=${OutDirNow//\//.}
JobFileName=${JobBaseName}.sh
JobFilePost=${JobBaseName}.2.sh

# Delete output file
[ -w $JobFileName    ] && rm $JobFileName
[ -w $JobFilePost    ] && rm $JobFilePost

#set -x           #Log all steps - debug
set -e            #Fail on any errors
shopt -s nullglob #Return nothing if no match

# Walk each of the ZV directories
for d in $InDirNow/*
do
    Dir=${d##*/}
    if [ "$Dir" != log ] && [ -d $d ]
    then
	Spec=$Dir
	SpecGrp=$Spec
	unset PGroup
	if [ "${Spec: -2}" == "p2" ]
	then
	    PGroup=${Spec: -2}
	    Spec=${Spec:0:-2}
	fi
	# Fit each file separately, with timeslices determined by deltaT
	unset NeedSummary
	# Covariance must be frozen. Ratios defined as ensemble average
	# ... so we don't have binned data to unfreeze covariance on each ensemble
	Cmd="MultiFit --freeze --delta 3 --iter 10000 --opnames"
	Cmd="$Cmd -i $InDir/$Dir/ -o $SpecGrp/"
	for p in $d/ZV_*_dt_*.h5
	do
	    FileName=${p##*/}
	    parts=(${FileName//_/ })
	    if (( ${#parts[@]} > 4 )) && [ "${parts[2]}" == "dt" ]
	    then
		DeltaT=${parts[3]}
		case "$DeltaT" in
		    32) FitRange=4:17:12:12;;
		    28) FitRange=4:15:10:10;;
		    24) FitRange=4:13:8:8;;
		    20) FitRange=4:11:6:6;;
		    16) FitRange=4:9:4:4;;
		    *)  unset FitRange; echo "Unrecognised DeltaT $DeltaT";;
		esac
		if [ -v FitRange ]
		then
		    echo "$Cmd -t $FitRange ${FileName},const" >> $JobFileName
		    NeedSummary=1
		fi
	    fi
	done
	# Now summarise the output
	if [ -v NeedSummary ]
	then
	    Cmd="FitSummary -i $SpecGrp/ -o $SpecGrp/Summary/ '*.model.*.h5'"
	    echo "$Cmd" >> $JobFilePost
	fi
    fi
done

# Job file needs to be executable
chmod u+x $JobFileName
