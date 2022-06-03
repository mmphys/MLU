#!/usr/bin/env bash
. common_utility.sh

if [[ "$pcorr" == "" ]]; then pcorr='../corr/2ptp2'; echo "default pcorr=$pcorr"; fi

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ "$mmplotfile_type" == "model" && "$mmplotfile_ext" == "h5" ]]; then #Silently skip pdfs

[[ -v srcmax ]] || srcmax=$((${#mmplotfile_ops[@]} - 1))
[[ -v srcmin ]] || srcmin=0
[[ -v snkmax ]] || snkmax=0
[[ -v snkmin ]] || snkmin=0
[[ -v ti || -v tf ]] && XRangeLabel="_${ti:-0}_${tf}"
[[ -v ti ]] && ti=$((ti-1)).8
[[ -v tf ]] && tf=${tf}.2

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark

FieldName="${field:-cosh}"
CorrPath="${pcorr}"
Basename="${mmplotfile_base}.${mmplotfile_corr_all}."
InBasename="${mmplotfile_path}".Basename
my_key="${key:-top right}"

#Function to return operator name
OpNameWords="${mmplotfile_ops[@]}"
OpName(n) = word( OpNameWords, n + 1 )

set term pdf
set key @my_key
set pointsize 0.6
set xlabel 'timeslice'
set xrange [${ti:-*}:${tf:-*}]
#set xrange [3.8:24.2]

set output Basename."${mmplotfile_ops_all}.${mmplotfile_type}_".FieldName."$XRangeLabel.${mmplotfile_seed}.pdf"
set title "Model vs data ${mmplotfile_corr}elated fit on [${mmplotfile_ti},${mmplotfile_tf}] ( ".OpNameWords." on 24^3 ensemble)"
set ylabel 'Synthesised correlator mass from fit'

set style fill transparent solid 0.2 noborder
plot \
  for [snk=${snkmin}:${snkmax}] for [src=${srcmin}:${srcmax}] \
    InBasename.OpName(snk)."_".OpName(src).".bootstrap.${mmplotfile_seed}.${mmplotvar_txt}" \
    using "t":(column(FieldName."_low")):(column(FieldName."_high")) with filledcurves title "Fit error ".OpName(snk)." - ".OpName(src), \
  for [snk=${snkmin}:${snkmax}] for [src=${srcmin}:${srcmax}] \
    InBasename.OpName(snk)."_".OpName(src).".bootstrap.${mmplotfile_seed}.${mmplotvar_txt}" \
    using "t":(column(FieldName)) with lines lw 0.5 title "Fit ".OpName(snk)." - ".OpName(src), \
  for [snk=${snkmin}:${snkmax}] for [src=${srcmin}:${srcmax}] \
    CorrPath."/${mmplotfile_base}_".OpName(snk)."_".OpName(src).".fold.${mmplotfile_seed}.${mmplotvar_txt}" \
  using "t":(column(FieldName)):(column(FieldName."_low")):(column(FieldName."_high")) with yerrorbars title "Data ".OpName(snk)." - ".OpName(src)
EOFMark

fi
done
