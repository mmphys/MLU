#!/bin/sh

if [[ "$pfold" == "" ]]; then pfold='../fold'; echo "default pfold=$pfold"; fi

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ "$mmplotfile_type" == "model" && "$mmplotfile_ext" == "h5" ]]; then #Silently skip pdfs

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark

FieldName="${field:-cosh}"
FoldPath="${pfold}"
Basename="${mmplotfile_base}.${mmplotfile_corr_all}."
InBasename="${mmplotfile_path}".Basename
my_key="${key:-top right}"

#Function to return operator name
OpName(n) = (n==0) ? "g5" : "gT5"
OpText(n) = (n==0) ? "Γ_5" : "Γ_4 Γ_5"

set term pdf
set key @my_key
set pointsize 0.6
set xlabel 'timeslice'
set xrange [${ti:-*}:${tf:-*}]
#set xrange [3.8:24.2]

set output Basename."${mmplotfile_ops_all}.${mmplotfile_type}_".FieldName.".${mmplotfile_seed}.pdf"
set title "Masses reconstructed from ${mmplotfile_corr}elated fit on [${mmplotfile_ti},${mmplotfile_tf}] ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)"
set ylabel 'Synthesised correlator mass from fit'

set style fill transparent solid 0.2 noborder
plot \
  for [snk=0:1] for [src=0:1] \
    InBasename.OpName(snk)."_".OpName(src).".bootstrap.${mmplotfile_seed}.${mmplotvar_dat}" \
    using "t":(column(FieldName."_low")):(column(FieldName."_high")) with filledcurves title "Error ".OpText(snk)." - ".OpText(src), \
  for [snk=0:1] for [src=0:1] \
    InBasename.OpName(snk)."_".OpName(src).".bootstrap.${mmplotfile_seed}.${mmplotvar_dat}" \
    using "t":(column(FieldName)) with lines lw 0.5 title "Midpoint ".OpText(snk)." - ".OpText(src), \
  for [snk=0:1] for [src=0:1] \
    FoldPath."/${mmplotfile_base}_".OpName(snk)."_".OpName(src).".fold.${mmplotfile_seed}.${mmplotvar_dat}" \
  using "t":(column(FieldName)):(column(FieldName."_low")):(column(FieldName."_high")) with yerrorbars title "Data ".OpText(snk)." - ".OpText(src)
EOFMark

fi
done
