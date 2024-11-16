#!/usr/bin/env bash

# Make Correlator links
# I.e. choose correlators where momentum applied to lighter quark
RootDir=analyse
RawCorr=corr_raw
Corr=corr

function GetWeight()
{
    case ${1:0:1} in
	h) Weight=3;;
	s) Weight=2;;
	l) Weight=1;;
	*) Weight=0;;
    esac
}

function IsMomentumZero()
{
    local p=$1
    local pstring=${p%%_*}
    p=${p#*_}
    Momentum=0
    if [ "$pstring" == p  ] && [ "${p:0:5}" != "0_0_0" ] ||
       [ "$pstring" == p2 ] && [ "${p%%_*}" != "0" ]
    then
	Momentum=1
    fi
}

if ! [ -d $RootDir/$RawCorr ]; then
    echo "Raw correlator directory missing: $RootDir/$RawCorr"
else
    cd $RootDir
    if [ -d $Corr ]; then
	echo "Deleting existing correlator links in $Corr"
	rm -rf $Corr
    fi
    mkdir $Corr
    cd $Corr
    RawCorr=../$RawCorr
    for d in $RawCorr/*; do
	s=${d##*/}
	if [ "${s:0:1}" != "2" ]; then
	    ln -s $d
	else
	    mkdir $s
	    cd $s
	    for f in ../$d/*; do
		# Get filename (without directory)
		g=${f##*/}
		# Extract quark names
		qtwo=${g%%_*}
		h=${g#*_}
		qone=${h%%_*}
		h=${h#*_}
		# See whether the momentum is zero
		IsMomentumZero $h
		# Get quark weights
		GetWeight $qtwo
		WeightTwo=$Weight
		GetWeight $qone
		# Create links for a momentum combination we like
		unset LinkForward
		if (( Momentum == 0 ))
		then
		    # Zero momentum
		    LinkForward=1
		    if (( WeightTwo != Weight )); then
			# Different weights -> link the reversed name to the same file
			ln -s $f ${qone}_${qtwo}_$h
		    fi
		else
		    # Non-zero momentum
		    if (( WeightTwo >= Weight )); then
			# Good: i.e. momentum was on the lightest quark
			LinkForward=1
		    else
			# Bad: i.e. momentum was on the heaviest quark
			# Link this name to the reversed file
			ln -s ${f%/*}/${qone}_${qtwo}_$h $g
		    fi
		fi
		if [ -v LinkForward ]; then
		    ln -s $f
		fi
	    done
	    cd ..
	fi
    done
fi
