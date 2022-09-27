#!/usr/bin/env bash

# Make wildcard failures return empty strings (not the unmodified wildcard string)
shopt -s nullglob
set -e
#set -x

# Folder containing the original bootstraps
Boot=bootstrap

# Make a file to submit Correlator folding / real / imaginary
BootRun=../$Boot/
BootNow=analyse/$Boot/
Cmd="corr -i $BootRun"

# 2pt functions
echo "${Cmd}2ptp2/'*'_ -o 2ptp2/ 'g5[PW]_g5[PW].*.h5',per"
echo "${Cmd}2ptp2/'*'_ -o 2ptp2/ 'gT5P_g5P.*.h5',nor"

# 3pt functions
for d in $BootNow/3*p2; do
    Spec=${d##*/}
    CmdBase="${Cmd}${Spec}/'*'_ -o ${Spec}/"
    echo "$CmdBase   gT_dt_'*'.h5,r"
    echo "$CmdBase gXYZ_dt_'*'.h5,i"
done

# Ward identities
echo "${Cmd}Ward/ -o Ward/ PJ5q_'*.h5',per"
echo "${Cmd}Ward/ -o Ward/ PA0_'*.h5',r"
