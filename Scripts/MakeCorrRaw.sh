#!/usr/bin/env bash

# Make a file to submit Correlator folding / real / imaginary
InBase=../bootstrap/
Cmd="corr -i $InBase"

# 2pt functions
echo "${Cmd}2ptp2/'*'_ -o 2ptp2/ 'g5[PW]_g5[PW].*.h5',per"

# 3pt functions
for Spec in 3pt_{h385,l,s}; do
    CmdBase="${Cmd}${Spec}p2/'*'_ -o ${Spec}p2/"
    echo "$CmdBase   gT_dt_'*'.h5,r"
    echo "$CmdBase gXYZ_dt_'*'.h5,i"
done

# Ward identities
echo "${Cmd}Ward/ -o Ward/ PJ5q_'*.h5',per"
echo "${Cmd}Ward/ -o Ward/ PA0_'*.h5',r"
