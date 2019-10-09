#!/bin/sh
function plot_one()
{
  local cmd="$1"
  if [[ $cmd != "" ]] ; then
    cmd="gnuplot -e \"${cmd}\""
    #echo $cmd
    eval "${cmd//$/\\$}"
  fi
}

function plot_file()
{
  local file=$1
  local file_array=(${file//./ })
  local corr=${file_array[0]}
  local corr_escaped=${corr//\_/\\\_}
  local type=${file_array[1]}
  local seed=${file_array[2]}
  local ext=${file_array[3]}
  local cmd="plot '${file}'"
  if [[ $ext == "txt" ]] ; then
    case $type in
      corr) local cmd3="set title 'Im( Correlator ) ${corr_escaped} (seed=${seed})'; ${cmd} using 1:4:5 with yerrorbars"
            local cmd2="set title 'Correlator ${corr_escaped} (seed=${seed})'; ${cmd} using 1:2:3 with yerrorbars"
            cmd="set title 'Correlator ${corr_escaped} (seed=${seed})'; set logscale y; ${cmd} using 1:(abs(\$2)):3 with yerrorbars";;
      mass | cosh | sinh | mass3pt ) cmd="set title '${type} ${corr_escaped} (seed=${seed})'; ${cmd} using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars";;
      *) echo "Unsupported file format: ${type}"
         cmd="set title '${type} ${corr_escaped} (seed=${seed})'; ${cmd}";;
    esac
  fi
  plot_one "$cmd"
  plot_one "$cmd2"
  plot_one "$cmd3"
}

while (( "$#" )); do
  echo About to plot $1
  plot_file $1
shift
done
