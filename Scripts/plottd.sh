#!/usr/bin/env bash
. common_utility.sh

###################################################
# Plot theory - data from fit
###################################################

PlotFunction()
{
gnuplot <<-EOFMark

#Command-line options
PlotFile="$PlotFile"
FieldNames="${fields:-log}"
FieldNamesText="${fieldtext:-effective mass}"
PDFSize="${size}"
Save="${save}"
xAxis="${x:-t}"

NumFields=words(FieldNames)

if( Save ne "" ) {
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  set output Save.'.pdf'
}

set key bottom right

fld=1
#set xrange [8.3:10.7]
#set xrange [7.8:11.2]
#set logscale y

plot \
  PlotFile using (stringcolumn('field') eq word(FieldNames,fld) ? column(xAxis) : NaN ) \
    :(column('theory_low')):(column('theory_high')) \
      with filledcurves notitle lc "skyblue" fs transparent solid 0.5,\
  '' using (stringcolumn('field') eq word(FieldNames,fld) ? column(xAxis) : NaN ):(column('theory')) \
      with linespoints notitle lc "blue" dt 5 pt 4,\
  '' using (stringcolumn('field') eq word(FieldNames,fld) ? column(xAxis) : NaN ) \
    :(column('data')):(column('data_low')):(column('data_high')) \
    with yerrorbars title FieldNamesText lc "red" pt 13

EOFMark
}

if [[ "$*" == "" ]];
then
  echo "$0"
  echo "Plot theory / data from fit."
  exit 2
fi

PlotFile="$*"
PlotFunction
