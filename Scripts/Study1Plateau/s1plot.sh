#!/bin/sh
FileParams=$1
if [[ "$FileParams" == "" ]]
then
  FileParams="$HOME/data/201910Plat/fit/Z2/h1_l_p_0_0_0.corr.g5_gT5.params.4147798751.txt"
  echo "FileParams defaulting to ${FileParams}"
fi

gnuplot -e "FileParams='${FileParams}'" ${0%.*}.gp
gnuplot -e "FileParams='${FileParams}'" ${0%.*}A.gp
