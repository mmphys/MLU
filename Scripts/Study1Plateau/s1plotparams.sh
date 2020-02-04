#!/bin/sh

# Loop through all the params files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ "$mmplotfile_type" != "params" ]]
then
  echo "Not params file: $mmplotfile_name"
else

if [[ "$chi" == "" ]]; then do_chi=0; else do_chi=1; fi

###################################################
# Plot the extracted energy levels and chi squared / d.o.f. associated with the model
###################################################

gnuplot <<-EOFMark

do_chi=${do_chi:-0}
chi_max=${chi:-0}
NumExp=${e:=1}

set term pdf

set pointsize 0.6
set xlabel 'initial fit time'
set xrange [${ti:=*}:${tf:=*}]
#set palette defined ( 15 "blue", 16 "red", 17 "green")

sChiDescr='{/Times:Italic χ}^2 per d.o.f.'
if( do_chi ) { sChiDescr=sChiDescr." <= ".sprintf("%g",chi_max) }
FitName=' from '.NumExp."-exponential ${mmplotfile_corr}elated fit using ${mmplotfile_ops_all//_/ }"
OutBase="${mmplotfile_base}.${mmplotfile_corr_all}.${mmplotfile_ops_all}."
OutSuffix=".${mmplotfile_seed}.pdf"

do for [MyFieldNum = 0:NumExp - 1] {
  MyField="E".MyFieldNum
  OutFile=OutBase.MyField.OutSuffix
  set output OutFile
  set label 1 OutFile noenhanced at screen 1, 0 right font "Arial,8" front textcolor "grey40" offset character -1, character 0.75
  if (MyFieldNum==0) {
    set arrow from graph 0, first 0.99656 to graph 1, first 0.99656 nohead front lc rgb "gray40" lw 0.25  dashtype "-"
    set label 2 "E_0=0.99656(95), ArXiv:1812.08791" at screen 0, 0 font "Arial,8" front textcolor "grey40" offset character 1, character 0.75
    #set yrange[0.95:1.05]
  } else {
    unset arrow
    unset label 2
    #unset yrange
    #set yrange[1.4:2.4]
  }
  MyTitle = MyField.FitName
  if( do_chi ) { MyTitle=MyTitle." (".sChiDescr.")" }
  set title MyTitle
  set ylabel MyField.'(t)'
  
  plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):(do_chi==0 ? column(MyField) : column("ChiSqPerDof") <= chi_max ? column(MyField) : 1/0 ):(column(MyField."_low")):(column(MyField."_high")) with yerrorbars title columnhead(1)
  #plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):(column(MyField)):(column(MyField."ErLow")):(column(MyField."ErHigh")) with yerrorbars title columnhead(1)
}

unset arrow
unset label 2
unset yrange
OutFile=OutBase.'chisq'.OutSuffix
set output OutFile
set label 1 OutFile noenhanced at screen 0.95, 0 right font "Arial,8" front textcolor "grey40" offset character 0, character 0.75
set title sChiDescr." ".FitName
unset xrange
#set xrange [3.8:12.2]
#set xrange [5.8:12.2]
set ylabel sChiDescr

MyField="ChiSqPerDof"
plot for [idx=0:*] "$PlotFile" index idx using (column("ti")-0.1+0.05*idx):(do_chi==0 ? column(MyField) : column("ChiSqPerDof") <= chi_max ? column(MyField) : 1/0 ):(column(-2)+1):(column(-2)+1) with linespoints pt variable lc variable title columnhead(1)
EOFMark

fi
done
