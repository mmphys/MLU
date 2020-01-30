#!/bin/sh

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile"
if [[ "$mmplotfile_ext" == "txt" ]]; then #Silently skip non-text files

###################################################
# Make a plot of the forward and backward-propagating waves
###################################################

gnuplot <<-EOFMark
#Command-line options
nt=${nt:=64}
ti=${ti:=1}
tf=${tf:=nt/2-1}
my_yrange="${yrange:=*:*}"
FieldName="${field:=y}"
RefVal=${ref:=-777}

#Full path to file
FileName="$PlotFile"
#filename (without path)
FileTitle="${mmplotfile_name}"
#File type (corr, cosh, mass, etc)
FileType="${mmplotfile_type}"
#Where to write the file
OutFile="${mmplotfile_base}.fb.${mmplotfile_type}.${mmplotfile_seed}.pdf"

set term pdf
set output OutFile
set xrange[ti:tf]
set yrange[@my_yrange]
set title "Forward & backward waves ".FileTitle noenhanced
set key top center

if( RefVal != -777 ) {
  set arrow from graph 0, first RefVal to graph 1, first RefVal nohead front lc rgb "gray40" lw 0.25  dashtype "-"
  set label "E_0=".sprintf("%f",RefVal) at graph 0, first RefVal font "Arial,8" front textcolor "grey40" offset character 1.5, character 0.35
}

ColName="(column(\"".FieldName."\"))"
ColLow="(".ColName." - column(\"".FieldName."_low\"))"
ColHigh="(".ColName." + column(\"".FieldName."_high\"))"

if( FileType eq "corr" ) {
  set logscale y
  ColName="(abs(".ColName."))"
}

plot FileName using "t":@ColName:@ColLow:@ColHigh with yerrorbars title "Forward", \
     '' using (nt - column("t")):@ColName:@ColLow:@ColHigh with yerrorbars title "Backward"

EOFMark

fi
done
