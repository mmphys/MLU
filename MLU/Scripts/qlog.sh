#!/usr/bin/env bash

#set -x
set -e

# For every entry in Qlog show the corresponding job log
for f in Qlog/*.out.log
do
    n="$(grep ' node jobs in ' "$f" || :)"
    if [ "$n" == "" ]; then
	n=0
    else
	n="${n#Looking for }"
	n=${n%% *}
    fi
    j=${f%.out.log}
    j=${j##*.}
    unset JobStat
    unset Msg
    if [ $(squeue -h -j $j 2> /dev/null|wc -l) == 1 ]; then
	JobStat="$(squeue -h -o %T -j $j)"
    fi
    d="$(grep 'Preparing to execute' "$f" || :)"
    if [ "$d" == "" ]
    then
	MsgPrefix="$f"
	d=$(tail -n 1 "$f")
	if [ "$d" == "No work to do" ]; then
	    JobStat="${JobStat:+${JobStat} }$d"
	else
	    JobStat="${JobStat:+${JobStat} }Error dispatching job"
	fi
    else
	d=${d##* }
	c=${d##*/}
	d=${d%/*}
	MsgPrefix="$d $c"
	g="$d/$c/log/$d.$c.$j.log"
	if ! [ -f "$g" ]; then
	    Msg="(log doesn't exist) $g"
	else
	    if grep 'Grid is finalizing now' $g &> /dev/null; then
		JobStat="${JobStat:+${JobStat} }Good"
	    fi
	    Msg="$(ls -laF "$g")"
	fi
    fi
    echo $n $j $MsgPrefix"${JobStat:+ (${JobStat})}""${Msg:+ $Msg}"
done
