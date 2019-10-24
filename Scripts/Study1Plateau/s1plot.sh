#!/bin/sh
PlotFile=$1
if [[ "$PlotFile" == "" ]]
then
  PlotFile="$HOME/data/201910Plat/fit/Z2/h1_l_p_0_0_0.corr.g5_gT5.params.4147798751.txt"
  #echo "FileParams defaulting to ${PlotFile}"
fi

# Get the directory name PlotFile is in
PlotDir=${PlotFile%/*}
if [[ "$PlotDir" == "" ]]
then
  if [[ "${PlotFile:0:1}" == "/" ]]
  then
    PlotDir=${PlotDir}/
  fi
elif [[ "${PlotDir: -1}" != "/" ]]
then
  PlotDir=${PlotDir}/
fi

###################################################
# Plot the extracted energy levels and chi squared / d.o.f.
###################################################

gnuplot <<-EOFMark
#show colornames
set term pdf

set pointsize 0.6
set xlabel 'initial fit time'
set xrange [4.8:11.2]
#set palette defined ( 15 "blue", 16 "red", 17 "green")

do for [MyFieldNum = 0:1] {
  MyField="E".MyFieldNum
  set output "a_fit_".MyField.".pdf"
  if (MyFieldNum==0) {
    set arrow from 4.8,0.99565 to 11.2,0.99565 nohead front lc rgb "gray40" lw 0.25  dashtype "-"
    set label "E_0=0.99656(95), ArXiv:1812.08791" at 4.9,0.99656 font "Arial,8" front textcolor "grey40" offset character 0,0.2
    TitleFieldName=MyField
  } else {
    unset arrow
    unset label
    TitleFieldName="Δ".MyField
  }
  set title TitleFieldName.' from 2-exponential fit ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
  set ylabel TitleFieldName.'(t)'
  plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):MyField:(column(MyField)-column(MyField."ErLow")):(column(MyField)+column(MyField."ErHigh")) with yerrorbars title columnhead(2)
}

set output "a_fit_chisq.pdf"
set title '{/Times:Italic χ}^2 per d.o.f. dependence on initial/final fit times ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
set xrange [4.8:12.2]
set ylabel 'Chi squared per degree of freedom'
plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):"ChiSqPerDof":(column(-2)+1):(column(-2)+1) with linespoints pt variable lc variable title columnhead(2)
EOFMark

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark
set term pdf

#Function to return operator name
OpName(n) = (n==0) ? "g5" : "gT5"
OpText(n) = (n==0) ? "Γ_5" : "Γ_4 Γ_5"
YOffset(snk,src) = (snk==0 && src==1) ? -.002 : 0

Seq="4147798751"
ti=6
tf=16
prefix="h1_l_p_0_0_0"
suffix=".txt"
InBasename=prefix.".corr_".ti."_".tf."."

set pointsize 0.6
set xlabel 'timeslice'
set xrange [4.8:11.2]

set output "a_reconstruct_".ti."_".tf.".pdf"
set title 'Masses reconstructed from fit on ['.ti.','.tf.'] ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
set xrange [4.8:12.2]
set ylabel 'Synthesised correlator mass from fit'

if (YOffset(0,1)!=0) {
  set label OpText(0)."-".OpText(1)." offset to avoid overlap with ".OpText(1)."-".OpText(0) \
    at graph 0,0 offset character 1,1.1 left font "Arial,8" textcolor "grey40"
}

set style fill transparent solid 0.2 noborder
plot for [snk=0:1] for [src=0:1] "$PlotDir".InBasename.OpName(snk)."_".OpName(src).".mass.".Seq.suffix \
  using 1:(\$2-\$3+YOffset(snk,src)):(\$2+\$4+YOffset(snk,src)) with filledcurves title "Error ".OpText(snk)." - ".OpText(src), \
  for [snk=0:1] for [src=0:1] "$PlotDir".InBasename.OpName(snk)."_".OpName(src).".mass.".Seq.suffix \
  using 1:(\$2+YOffset(snk,src)) with lines lw 0.5 title "Midpoint ".OpText(snk)." - ".OpText(src)
EOFMark
