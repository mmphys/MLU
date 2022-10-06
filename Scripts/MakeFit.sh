#!/bin/env bash

shopt -s nullglob # globs return empty strings if no match

# Make an analysis job which will perform fits
# Optional environment variables:
#   uncorr= Perform uncorrelated fits
#   minuit= Use Minuit2 as fitter
#   num     Number of bootstrap replicas to fit (default all)
JobRoot=${JobRoot:-analyse}
JobIn=${JobIn:-corr}
JobOut=${JobOut:-fit}

# How to reach correlators during run
CorrDir=../$JobIn
# How to reach directories now (while this script runs)
CorrDirNow=$JobRoot/$JobIn
FitDirNow=$JobRoot/$JobOut

# Don't expect these to be overriden
JobBaseName=${FitDirNow//\//.}
JobFileName=${JobBaseName}.sh
JobFilePost=${JobBaseName}.2.sh
Point=g5P_g5P
Wall=g5P_g5W

if ! [ -r ./FitRanges.sh ]
then
    echo "FitRanges.sh doesn't exist"
    exit 1
fi
. ./FitRanges.sh
if ! [[ $(type -t FitRanges) == function ]]
then
    echo "FitRanges.sh did not define function FitRanges()"
    exit 1
fi

# Get a relative weight for a quark
function GetWeight()
{
    case ${1:0:1} in
        h) Weight=3;;
        s) Weight=2;;
        l) Weight=1;;
        *) Weight=0;;
    esac
}

# Get the subset of q2_q1 combinations we're interested in
function GetCombos()
{
    local Dir=$1
    local qone
    local qtwo
    local Weight
    local WeightTwo
    for f in $Dir/*.h5; do
	# Grab quark names from first two fields
	g=${f##*/}
	qtwo=${g%%_*}
	g=${g#*_}
	qone=${g%%_*}
	# Get quark weights
	GetWeight $qtwo
	WeightTwo=$Weight
	GetWeight $qone
	# Ignore Q1 heavier than Q2
	if (( Weight <= WeightTwo )); then
	    Combos[${qtwo}_$qone]=${qtwo:0:1}${qone:0:1}
	fi
    done
}

# Delete output files
[ -w $JobFileName    ] && rm $JobFileName
[ -w $JobFilePost    ] && rm $JobFilePost

for Dir in 2ptp2
do
    # Directories relative to when job runs
    unset Combos
    declare -A Combos
    InDir=$CorrDir/$Dir
    InDirNow=$CorrDirNow/$Dir
    GetCombos $InDirNow
    echo "$Dir Combos: ${!Combos[@]}"
    for Combo in ${!Combos[@]}
    do
	Out=$Dir/$Combo
	unset Range1
        unset Range2
        unset Options
	FitRanges ${Combos[$Combo]}
	if ! [ -v Range1 ] || ! [ -v Range2 ]
        then
            echo "Unsupported job type ${Combos[$Combo]}"
            exit 1
        fi
	for p in $InDirNow/${Combo}_*_${Point}*h5
	do
	    Prefix=${p##*/}               #Leave the filename only
	    Prefix=${Prefix%%${Point}*h5} #Chop off everything past point/wall
	    Cmd="MultiFit --iter 10000 -e 2 --mindp 8"
	    if [ -v uncorr  ]; then Cmd="$Cmd --uncorr"; fi
	    if [ -v minuit  ]; then Cmd="$Cmd --fitter minuit2"; fi
	    if [ -v num     ]; then Cmd="$Cmd -n $num"; fi
	    if [ -v Options ]; then Cmd="$Cmd $Options"; fi
	    # NB the trailing 'n' = normalise by energy
	    Cmd="$Cmd -o $Out/ -i $InDir/"
	    #Cmd="$Cmd ${Prefix}{${Point},${Wall}}*h5,,,n"
	    Cmd="$Cmd ${Prefix}${Point}*h5,t=${Range1},eNorm=true"
	    Cmd="$Cmd ${Prefix}${Wall}*h5,t=${Range2},eNorm=true"
	    echo "$Cmd" >> $JobFileName
	done
	# Now summarise the output
	Cmd="FitSummary -o $Dir/Summary/ -i $Out/ ${Combo}_'*.model.*'.h5"
	echo "$Cmd" >> $JobFilePost
    done
done

# Job file needs to be executable
chmod u+x $JobFileName
