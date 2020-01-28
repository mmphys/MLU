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
  local file_short=$file
  local ext=${file_short//*.}
        file_short=${file_short%.*}
  local seed=${file_short//*.}
        file_short=${file_short%.*}
  local type=${file_short//*.}
  local corr=${file_short%.*}
  local corr_escaped=${corr//\_/\\\_}
  local cmd="plot '${file}'"
# cmd="set xrange [5.9:16.1]; ${cmd}"
  cmd="set xrange [2:62]; ${cmd}"
#  cmd="set term pdf; set output '${file%.*}.pdf'; ${cmd}"
  if [[ $ext == "txt" ]] ; then
    local title=" ${corr_escaped} (seed=${seed})"
    local fields="2:(\$2-\$3):(\$2+\$4) with yerrorbars"
    case $type in
      corr) title="set title 'Correlator ${title}'"
            plot_one "${title}; ${cmd} using 1:6:(\$6-\$7):(\$6+\$8) with yerrorbars title 'imaginary' lc rgb 'red', '${file}' using 1:${fields} title 'real' lc rgb 'blue'"
            cmd="${title}; set logscale y; ${cmd} using 1:(abs(\$2)):(abs(\$2)-\$3):(abs(\$2)+\$4) with yerrorbars title 'real' lc rgb 'blue'";;
      mass | cosh | sinh ) cmd="set title '${type} ${title}'; ${cmd} using 1:${fields} title '${type}' lc rgb 'red'";;
      *) echo "Unsupported file format: ${type}"
         cmd="set title '${type} ${title}'; ${cmd}";;
    esac
  fi
  plot_one "$cmd"
}

while (( "$#" )); do
  #echo About to plot $1
  plot_file $1
shift
done
