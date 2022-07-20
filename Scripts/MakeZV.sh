#!/bin/env bash

# Make an analysis job which will perform fits
# Optional environment variables:
JobRoot=${JobRoot:-analyse}
InCorr=${InCorr:-corr}
InFit=${InFit:-fit}
Out=${Out:-ZV}

# How to reach correlators during run
CorrDir=../$InCorr
FitDir=../$InFit
# How to reach directories now (while this script runs)
CorrDirNow=$JobRoot/$InCorr
FitDirNow=$JobRoot/$InFit
OutDirNow=$JobRoot/$Out

# Don't expect these to be overriden
JobBaseName=${OutDirNow//\//.}
JobFileName=${JobBaseName}.sh

# Delete output file
[ -w $JobFileName    ] && rm $JobFileName

#set -x           #Log all steps - debug
set -e            #Fail on any errors
shopt -s nullglob #Return nothing if no match

# Walk each of the 3pt directories
for d in $CorrDirNow/3pt_*
do
    Dir=${d##*/}
    Spec=${Dir:4}
    SpecGrp=$Spec
    unset PGroup
    if [ "${Spec: -2}" == "p2" ]
    then
	PGroup=${Spec: -2}
	Spec=${Spec:0:-2}
    fi
    # Fit lists must have been chosen
    FitList="Fit_${SpecGrp}.txt"
    if ! [ -r "$FitDirNow/$FitList" ]
    then
	echo "$SpecGrp: $FitDirNow/$FitList doesn't exist"
    else
	echo "$SpecGrp: $FitDirNow/$FitList exists"
	# Get a list (associative array) of all the quarks
	unset Q
	declare -A Q
	for f in $d/*_gT_dt_*_p2_0_ps2_0_*.h5
	do
	    g=${f##*/}
	    parts=(${g//_/ })
	    if (( ${#parts[@]} > 3 )) && [ "${parts[1]}" == "${parts[2]}" ]
	    then
		Q[${parts[1]}]=${parts[0]}
	    fi
	done
	# Now make a job for each Quark
	Cmd="CRatio --i3 $CorrDir/$Dir/ --i2 $CorrDir/2pt${PGroup}/ -o $SpecGrp/"
	Cmd="$Cmd --efit $FitDir/Fit_${SpecGrp}.txt"
	unset f
	for q in ${!Q[@]}
	do
	    echo "$Cmd '*'_${q}_${q}_gT_dt_'*'_p2_0_ps2_0_'*'.h5" >> $JobFileName
	done
    fi
done

# Job file needs to be executable
chmod u+x $JobFileName
