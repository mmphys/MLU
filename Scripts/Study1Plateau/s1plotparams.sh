#!/bin/sh

# Loop through all the params files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile"
if [[ "$mmplotfile_type" != "params" ]]
then
  echo "Not params file: $mmplotfile_name"
else

###################################################
# Plot the extracted energy levels and chi squared / d.o.f. associated with the model
###################################################

gnuplot <<-EOFMark
#show colornames
set term pdf

set pointsize 0.6
set xlabel 'initial fit time'
#  set xrange [4.8:9.2]
#set palette defined ( 15 "blue", 16 "red", 17 "green")

do for [MyFieldNum = 0:1] {
  MyField="E".MyFieldNum
  set output "${mmplotfile_base}.${mmplotfile_corr}.".MyField.".${mmplotfile_seed}.pdf"
  if (MyFieldNum==0) {
    set arrow from 4.8,0.99656 to 11.2,0.99656 nohead front lc rgb "gray40" lw 0.25  dashtype "-"
    set label "E_0=0.99656(95), ArXiv:1812.08791" at 4.9,0.99656 font "Arial,8" front textcolor "grey40" offset character 0,0.2
    TitleFieldName=MyField
  } else {
    unset arrow
    unset label
    set yrange[0.9:1.9]
    TitleFieldName="Δ".MyField
  }
  set title TitleFieldName.' from 2-exponential ${mmplotfile_corr}elated fit ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
  set ylabel TitleFieldName.'(t)'
  plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):MyField:(column(MyField)-column(MyField."ErLow")):(column(MyField)+column(MyField."ErHigh")) with yerrorbars title columnhead(1)
}

unset yrange
set output "${mmplotfile_base}.${mmplotfile_corr}.chisq.${mmplotfile_seed}.pdf"
set title '{/Times:Italic χ}^2 per d.o.f. dependence on initial/final ${mmplotfile_corr}elated fit times ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
#  set xrange [3.8:12.2]
#  set xrange [5.8:12.2]
set ylabel 'Chi squared per degree of freedom'
plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):"ChiSqPerDof":(column(-2)+1):(column(-2)+1) with linespoints pt variable lc variable title columnhead(1)
EOFMark

fi
done
