#!/usr/bin/env bash

# Command-line arguments are a list of directories to show run summary for
# Optional environment variables:
# grepfor: grep command e.g. "-B 2 Converged", default 'Grid is finalizing now'
# lines:   number of lines at end of file to display

#set -x
set -e
shopt -s nullglob
#shopt -s globstar

function ShowMe()
{
    for f in "$@"; do echo $f; for c in $f/*; do
	n="$(ls -1R $c/ | wc -l)"
	for g in $c/log/*[0-9].[0-9][0-9][0-9][0-9][0-9].log; do
	    j=${g%.log}
	    c=${j%.*}
	    j=${j##*.}
	    c=${c##*.}
	    if [ "${lines+y}" == "y" ]; then
		echo "$c $j $n $(tail -n $lines $g)"
	    elif [ "${grepfor+y}" == "y" ]; then
		echo "$c $j $n $(grep $grepfor $g)"
	    else
		if [ $(squeue -h -j $j 2> /dev/null|wc -l) == 1 ]; then
		    JobStat="$(squeue -h -o %T -j $j)"
		elif grep 'Grid is finalizing now' $g &> /dev/null; then
		    JobStat=Good
		else
		    JobStat=Bad
		fi
		echo $c $j $n $JobStat $g
	    fi
	done
    done; done
}

if [ ${#@} == 0 ]; then
    ShowMe gf*
else
    ShowMe "$@"
fi
