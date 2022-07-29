#!/usr/bin/env bash

#set -x           #Log all steps - debug
set -e            #Fail on any errors
shopt -s nullglob #Return nothing if no match

# Expects to be run from ZVFit directory
# Makes ZV_Fit*.sh files choosing best fits by p-value
# Not all fits might be available, so edit the fit files once generated

function MakeAll()
{
    local Count=0
    for d in *
    do
	Dir=${d##*/}
	if [ "$Dir" != "log" ] && [ -d $d ]
	then
	    #Directory name is spectator and momentum group
	    Spec=$Dir
	    SpecGrp=$Spec
	    unset PGroup
	    if [ "${Spec: -2}" == "p2" ]
	    then
		PGroup=${Spec: -2}
		Spec=${Spec:0:-2}
	    fi
	    # Make sure output doesn't already exist
	    Outfile=ZVFit_${SpecGrp}.auto.txt
	    if [ -f $Outfile ]
	    then
		echo $Outfile already exists. Skipping.
	    elif [ -d $d/Summary ]
	    then
		echo "# ZV fits automatically selected" >> $Outfile
		echo "# $0" >> $Outfile
		echo "# $(date)" >> $Outfile
		echo "# Choose fit with most data points then best p-value" >> $Outfile
		# Now walk summary directory
		for f in $d/Summary/*.params_sort.*.txt
		do
		    Filename=${f##*/}
		    # sort file by NumDataPoints (field 5) reversed
		    #         then pValue sort sequence (field 1)
		    # grab ti and tf from the first (best) result
		    Fields=($(sort -V -k 5r -k 1 $f | awk '/^[-0-9]/ {print $2, $3;exit}'))
		    if [ "${#Fields[@]}" != 2 ]
		    then
			echo Bad $f
		    else
			IFS='.' read -ra Parts <<< "$Filename"
			Base=(${Parts[0]//_/ })
			Quark=${Base[1]}
			DeltaT=${Base[3]}
			Parts[1]=${Parts[1]}_${Fields[0]}_${Fields[1]}
			Parts[3]=model
			Parts[5]=h5
			Newname="${Parts[@]}"
			Newname=${Newname// /.}
			echo $Quark $DeltaT ${Parts[2]//_/ } $d/$Newname >> $Outfile
		    fi
		done
	    fi
	fi
    done
}

MakeAll
unset -f MakeAll
