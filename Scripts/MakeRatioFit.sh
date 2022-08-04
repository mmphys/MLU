#!/bin/env bash

# Make an analysis job which will perform fits of the ratios
# Optional environment variables:
JobRoot=${JobRoot:-analyse}
JobIn=${In:-ratio}
JobOut=${Out:-RatioFit}

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

# Walk each of the ratio directories
CmdSuffix=_g5P_g5P.fold.1835672416.h5,const
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
        unset NeedSummary
        Cmd="MultiFit --freeze --delta 1 --mindof 0 --iter 10000 --opnames"
        Cmd="$Cmd -i $InDir/$Dir/ -o $SpecGrp/ -t"
	# Do we have fit ranges
	FitRanges=FitRangesRatio_$SpecGrp.txt
	if ! [ -r $FitRanges ]
	then
	    echo "$FitRanges missing"
	else
	    echo "$FitRanges"
	    # Covariance must be frozen. Ratios defined as ensemble average
	    # ... so we don't have binned data to unfreeze covariance
	    while read ThisLine || [ -n "$ThisLine" ]
	    do
		# Echo comments to screen
		if [ "${ThisLine:0:1}" == "#" ]
		then
		    echo $ThisLine
		elif [ -n "$ThisLine" ]
		then
		    echo "$Cmd $ThisLine$CmdSuffix" >> $JobFileName
                    NeedSummary=1
		fi
	    done < $FitRanges
	fi
        # Now grab the midpoints of R1 and R2
        if [ -v dor12 ]
        then
            for f in $d/R[12]_*_g5P_g5P.*.h5
            do
                Filename=${f##*/}
                Parts=(${Filename//_/ })
                if (( ${#Parts[@]} == 10 ))
                then
                    DTHalf=$(( ${Parts[5]} / 2 ))
                    echo "$Cmd ${DTHalf}:${DTHalf} ${Filename},const" >> $JobFileName
                    NeedSummary=1
                fi
            done
        fi
        # Now summarise the output
        if [ -v NeedSummary ]
        then
            Cmd="FitSummary -i $SpecGrp/ -o $SpecGrp/Summary/"
            Cmd="$Cmd '*.model.*.h5'"
            echo "$Cmd" >> $JobFilePost
        fi
    fi
done

# Job file needs to be executable
chmod u+x $JobFileName
