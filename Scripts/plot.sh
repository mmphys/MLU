#!/bin/sh
file=$1

function do_plot()
{
  local cmd="$1"
  if [[ $cmd != "" ]] ; then
    cmd="gnuplot -e \"${cmd}\""
    #echo $cmd
    eval "${cmd//$/\\$}"
  fi
}

file_array=(${file//./ })
corr=${file_array[0]}
corr_escaped=${corr//\_/\\\_}
type=${file_array[1]}
seed=${file_array[2]}
ext=${file_array[3]}
cmd="plot '${file}'"
if [[ $ext == "txt" ]] ; then
  case $type in
corr) cmd2="set title 'Im( Correlator ) ${corr_escaped} (seed=${seed})'; ${cmd} using 1:4:5 with yerrorbars"
          cmd="set title 'Correlator ${corr_escaped} (seed=${seed})'; set logscale y; ${cmd} using 1:2:3 with yerrorbars";;
    mass | mass3pt) cmd="set title '${type} ${corr_escaped} (seed=${seed})'; ${cmd} using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars";;
    *) echo "Unsupported file format: ${type}"
       cmd="set title '${type} ${corr_escaped} (seed=${seed})'; ${cmd}";;
  esac
fi
do_plot "$cmd"
do_plot "$cmd2"
