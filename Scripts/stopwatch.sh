#!/usr/bin/env bash

formatSeconds()
{
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02d:%02d:%02d" $h $m $s
}

startTime=$(date +%s)

while [ 1 ]
do
    currentTime=$(date +%s)
    timePassed=$[$currentTime-startTime]

    echo -n "$(formatSeconds $timePassed)"

    sleep 0.5
    echo -en "\b\b\b\b\b\b\b\b"
done
