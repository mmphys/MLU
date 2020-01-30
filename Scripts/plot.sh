#!/bin/sh

if [[ "$pdf" == "1" ]]; then SaveFile=1; else SaveFile=0; fi

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile"
if [[ "$mmplotfile_ext" == "txt" ]]; then #Silently skip non-text files
log_limit=1
if [[ "$mmplotfile_type" == "corr" ]]; then log_limit=2; fi
for(( do_log=0 ; do_log < log_limit ; do_log = do_log + 1 ))
do
###################################################
# Make a plot of the forward and backward-propagating waves
###################################################

gnuplot <<-EOFMark
#Command-line options
my_xrange="${ti:=*}:${tf:=*}"
my_yrange="${yrange:=*:*}"
FieldName="${field:=y}"
RefVal=${ref:=-777}
do_log=${do_log:=0}
SaveFile=${SaveFile:=0}

#Full path to file
PlotFile="$PlotFile"
#Parsed bits of the filename
FileName="${mmplotfile_name}"
FileBase="${mmplotfile_base}"
FileType="${mmplotfile_type}"
FileSeed="${mmplotfile_seed}"
FileExt="${mmplotfile_ext}"
#Where to write the file
OutFile=FileBase.".".FileType
MyTitle=FileBase." (seed=".FileSeed.")"
if( do_log ) {
  OutFile=OutFile."_log"
  MyTitle=MyTitle." log"
}
OutFile=OutFile.".".FileSeed.".pdf"

if( SaveFile == 1 ) {
  set term pdf
  set output OutFile
  set pointsize 0.5
}

set xrange[@my_xrange]
set yrange[@my_yrange]

if( RefVal != -777 ) {
  set arrow from graph 0, first RefVal to graph 1, first RefVal nohead front lc rgb "gray40" lw 0.25  dashtype "-"
  set label sprintf("%f",RefVal) at graph 0, first RefVal font "Arial,8" front textcolor "grey40" offset character 1.5, character 0.35
}

#FieldNameY="column(\"".FieldName."\")"
#FieldNameLow="column(\"".FieldName."_low\")"
#FieldNameHigh="column(\"".FieldName."_high\")"
FieldNameLow=FieldName."_low"
FieldNameHigh=FieldName."_high"

AbsMin(y,low,high)=sgn(y) < 0 ? -(y+high) : y - low
AbsMax(y,low,high)=sgn(y) < 0 ? -(y-low) : y + high

if( FileType eq "corr" ) {
  # Correlator: Plot real and imaginary on non-log scale
  if( do_log ) { set logscale y }
  set title "Correlator ".MyTitle noenhanced
  plot PlotFile using 1:6:(\$6-\$7):(\$6+\$8) with yerrorbars title 'imaginary' lc rgb 'red', \
      '' using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars title 'real' lc rgb 'blue'
} else {
  #Not a correlator
  if( FileType eq "mass" ) {
    set title "Exponential mass ".MyTitle noenhanced
    plot PlotFile using "t":(abs(column(FieldName))):(AbsMin(column(FieldName),column(FieldNameLow),column(FieldNameHigh))):(AbsMax(column(FieldName),column(FieldNameLow),column(FieldNameHigh))) \
      with yerrorbars lc rgb 'blue' notitle
  } else {
    if( FileType eq "cosh" ) {
      MyTitle="Cosh mass ".MyTitle
    } else { MyTitle=FileType." ".MyTitle }
    set title MyTitle noenhanced
    plot PlotFile using "t":FieldName:(column(FieldName)-column(FieldNameLow)):(column(FieldName)+column(FieldNameHigh)) with yerrorbars lc rgb 'blue' notitle
  }
}

EOFMark
done; fi; done
