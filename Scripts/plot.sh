#!/bin/sh

PlotFunction()
{
# Loop through all the files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile"
if [[ "$mmplotfile_ext" == "txt" ]]; then #Silently skip non-text files
#log_limit=1
#if [[ "$mmplotfile_type" == "corr" ]]; then log_limit=2; fi
#for(( do_log=0 ; do_log < log_limit ; do_log = do_log + 1 ))
#do
###################################################
# Make a plot of the forward and backward-propagating waves
###################################################

gnuplot <<-EOFMark
#Command-line options
my_xrange="${ti:-*}:${tf:-*}"
my_yrange="${yrange:-*:*}"
my_key="${key:-top right}"
FieldNames="${fields:-cosh}"
RefVal=${ref:--777}
do_log=${log:-0}
SaveFile=${save:-0}
nt=${nt:-0}

#Full path to file
PlotFile="$PlotFile"
#Parsed bits of the filename
FileName_no_ext="${mmplotfile_name_no_ext}"
FileName="${mmplotfile_name}"
FileBase="${mmplotfile_base}"
FileType="${mmplotfile_type}"
FileSeed="${mmplotfile_seed}"
FileExt="${mmplotfile_ext}"

#Where to write the file (if specified)
OutFile=FileBase.".".FileType
MyTitle=FileName_no_ext
f = 1
while (f <= words(FieldNames) ) {
  OutFile=OutFile."_".word(FieldNames,f)
  f=f+1
}
fb_min=0
fb_max=0
array fb_prefix[2]
fb_prefix[1]=""
if ( nt != 0) {
  # We want the backward propagating wave
  fb_max = 1
  fb_prefix[2]="back "
  if ( nt < 0 ) {
    # We ONLY want the backward propagating wave
    nt = -nt
    fb_min = 1
    OutFile=OutFile."_b"
    MyTitle=MyTitle." backward"
  } else {
    fb_prefix[1]="fwd "
    OutFile=OutFile."_fb"
    MyTitle=MyTitle." forward/backward"
  }
}
if( do_log ) {
  OutFile=OutFile."_log"
  MyTitle=MyTitle." log"
}
OutFile=OutFile.".".FileSeed.".pdf"
#print OutFile

if( SaveFile == 1 ) {
  set term pdf
  set output OutFile
  set pointsize 0.5
}

set key @my_key
set xrange[@my_xrange]
set yrange[@my_yrange]
set title MyTitle noenhanced

if( RefVal != -777 ) {
  set arrow from graph 0, first RefVal to graph 1, first RefVal nohead front lc rgb "gray40" lw 0.25  dashtype "-"
  set label sprintf("%f",RefVal) at graph 0, first RefVal font "Arial,8" front textcolor "grey40" offset character 1.5, character 0.35
}

#FieldNameY="column(\"".FieldName."\")"
#FieldNameLow="column(\"".FieldName."_low\")"
#FieldNameHigh="column(\"".FieldName."_high\")"
#FieldNameLow=FieldName."_low"
#FieldNameHigh=FieldName."_high"

AbsMin(y,low,high)=sgn(y) < 0 ? -(high) : low
AbsMax(y,low,high)=sgn(y) < 0 ? -(low) : high

set linetype 1 lc rgb 'blue'
set linetype 2 lc rgb 'red'

if( do_log ) { set logscale y }

sForBack="Fwd Back"
plot for [fld in FieldNames] for [f=fb_min:fb_max] \
    PlotFile using (column(1) == 0 ? 0 : f==0 ? column(1) : nt - column(1)):(column(fld)):(column(fld."_low")):(column(fld."_high")) with yerrorbars title fb_prefix[f+1].fld

#if( FileType eq "corr" ) {
  # Correlator: Plot real and imaginary on non-log scale
#  set title "Correlator ".MyTitle noenhanced
#  plot PlotFile using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars title 'real' lc rgb 'blue', \
      '' using 1:6:(\$6-\$7):(\$6+\$8) with yerrorbars title 'imaginary' lc rgb 'red'
#} else {
  #Not a correlator
#  if( FileType eq "mass" ) {
#    set title "Exponential mass ".MyTitle noenhanced
#    plot PlotFile using "t":(abs(column(FieldName))):(AbsMin(column(FieldName),column(FieldNameLow),column(FieldNameHigh))):(AbsMax(column(FieldName),column(FieldNameLow),column(FieldNameHigh))) \
      with yerrorbars lc rgb 'blue' notitle
#  } else {
#    if( FileType eq "cosh" ) {
#      MyTitle="Cosh mass ".MyTitle
#    } else { MyTitle=FileType." ".MyTitle }
#    set title MyTitle noenhanced
#    plot PlotFile using "t":FieldName:(column(FieldName)-column(FieldNameLow)):(column(FieldName)+column(FieldNameHigh)) with yerrorbars lc rgb 'blue' notitle
#  }
#}

EOFMark
#done
fi; done
}

if [[ "$*" == "" ]];
then
  echo "$0"
  echo "Plot summary data."
  echo "Precede with optional modifiers (key=value):"
  echo "ti     Initial timeslice"
  echo "tf     Final timeslice"
  echo "key    location for key (default: top right)"
  echo "fields Names of fields to display (default: cosh)"
  echo "ref    y-value for reference line"
  echo "log    1 to plot y on log scale"
  echo "save   Save the plot to pdf (NB: filename will be auto generated)"
  echo "nt     Plot backward propagating wave as well (or backward only if nt<0)"
  exit 2
fi

if [[ "$pdf" == "1" ]]; then SaveFile=1; else SaveFile=0; fi
if [[ "$nt" != "" && "$tf" == "" ]]; then tf=$(( (nt < 0 ? -nt : nt)/2 )); fi

PlotFunction $*
