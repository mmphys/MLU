#!/usr/bin/env bash

#set -x
#set -e
shopt -s globstar

function ShowMe()
{
    # I want to sort all the log files from all directories in time order
    local Dirs="$@"
    if [ ${#@} != 1 ]; then Dirs="{${Dirs// /,}}"; fi
    Dirs="ls -Rt $Dirs/**/log/*[0-9].[0-9][0-9][0-9][0-9][0-9].log"
    #echo "\"$Dirs\""
    Dirs="$(shopt -s globstar; eval $Dirs)"
    #echo "\"$Dirs\""
    for g in $Dirs; do
	j=${g%.log}
	c=${j%.*}
	j=${j##*.}
	c=${c##*.}
	echo -e "====================\n$c $j $g"
	ls -laF $g
	grep sourceme ${g%.*}.sh
	grep 'mpirun -np' $g
	grep -A 2 '<BatchSize>' ${g%log}xml
	grep MLU $g
	tail -n 1 $g
	ls -laF Qlog/*$j*
	for l in Qlog/*$j*; do
	    echo -e "----------\n$l"
	    cat $l
	done
    done
}

# Command-line arguments are a list of directories to show run summary for
if [ ${#@} == 0 ]; then
    ShowMe semilep.run*
else
    ShowMe "$@"
fi
