#!/bin/sh

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile"
if [[ "$mmplotfile_type" == "model" && "$mmplotfile_ext" != "pdf" ]]; then #Silently skip pdfs

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark
set term pdf

#Function to return operator name
OpName(n) = (n==0) ? "g5" : "gT5"
OpText(n) = (n==0) ? "Γ_5" : "Γ_4 Γ_5"
YOffset(snk,src) = (snk==0 && src==1) ? -.002 : 0

set pointsize 0.6
set xlabel 'timeslice'
set xrange [${ti:=${mmplotdefault_ti}}-0.2:${tf:=${mmplotdefault_tf}}+0.2]
#set xrange [3.8:24.2]

Basename="${mmplotfile_base}.${mmplotfile_corr_all}."
InBasename="${mmplotfile_path}".Basename
set output Basename."${mmplotfile_ops_all}.${mmplotfile_type}.${mmplotfile_seed}.pdf"
set title "Masses reconstructed from ${mmplotfile_corr}elated fit on [${mmplotfile_ti},${mmplotfile_tf}] ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
set ylabel 'Synthesised correlator mass from fit'

if (YOffset(0,1)!=0) {
  set label OpText(0)."-".OpText(1)." offset to avoid overlap with ".OpText(1)."-".OpText(0) \
    at graph 0,0 offset character 1,1.1 left font "Arial,8" textcolor "grey40"
}

set style fill transparent solid 0.2 noborder
plot for [snk=0:1] for [src=0:1] InBasename.OpName(snk)."_".OpName(src).".mass.${mmplotfile_seed}.${mmplotvar_dat}" \
  using 1:(\$2-\$3+YOffset(snk,src)):(\$2+\$4+YOffset(snk,src)) with filledcurves title "Error ".OpText(snk)." - ".OpText(src), \
  for [snk=0:1] for [src=0:1] InBasename.OpName(snk)."_".OpName(src).".mass.${mmplotfile_seed}.${mmplotvar_dat}" \
  using 1:(\$2+YOffset(snk,src)) with lines lw 0.5 title "Midpoint ".OpText(snk)." - ".OpText(src), \
  for [snk=0:1] for [src=0:1] "${mmplotpath_corr}${mmplotfile_base}_".OpName(snk)."_".OpName(src).".mass.${mmplotfile_seed}.${mmplotvar_dat}" \
  using 1:(\$2+YOffset(snk,src)):(\$2+YOffset(snk,src)-\$3):(\$2+YOffset(snk,src)+\$4) with yerrorbars title "Data ".OpText(snk)." - ".OpText(src)
EOFMark

fi
done
