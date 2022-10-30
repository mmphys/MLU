#!/usr/bin/env bash

#set -x
set -e

export qSrc=h385
export qSnk=l
export qSpec=s

#pSnk=0 TI='11 12 14 12' TF='15 19 21 26' yrangeMEL=0.52:0.537 yrangeR3=0.75:1.05 FitMEL.sh
pSnk=1 TI='10 12 14 10' TF='13 16 20 26' yrangeMEL=0.48:0.496 yrangeR3=0.7:0.96 FitMEL.sh
