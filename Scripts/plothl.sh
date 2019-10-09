#!/bin/sh

if [[ "$1" == "" ]]
then
  inP2="0"
  echo "P2 defaulting to ${inP2}"
else
  inP2="$1"
fi

case $2 in
  0 | i) inVec="$2" ;;
  *) inVec=0; echo "Vector type defaulting to V${inVec}";;
esac

if [[ "$3" == "" ]]
then
  inDeltaT="16"
  echo "DeltaT defaulting to ${inDeltaT}"
else
  inDeltaT="$3"
fi

prefix="V${inVec}_${inDeltaT}_p${inP2}"
echo "Plotting ${prefix}"

gnuplot -e "prefix='${prefix}'; seq='4147798751'; qh='h1'; ql='l'; inVec='${inVec}'; inDeltaT='${inDeltaT}'; inP2='${inP2}'" $HOME/src/Utility/Scripts/plothl.gp
