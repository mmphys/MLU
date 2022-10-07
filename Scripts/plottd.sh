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
my_yrange="${yrange:-*:*}"
FieldNames="${fields:-log}"
FieldNamesText="${fieldtext:-effective mass}"
MyTitle="${title}"
PDFSize="${size:-6in,3in}"
Save="${save}"
xAxis="${x:-t}"
ModelMin=${mmin:-0}
ModelMax=${mmax:--1}

##################################################
# TODO: Convert to command-line parameters

# These are the files data points outside the fit range come from
OriginalNames='../F1M/corr/3pt_sp2/quark_l_h385_gT_dt_16_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5 ../F1M/corr/3pt_sp2/quark_l_h385_gT_dt_20_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5 ../F1M/corr/3pt_sp2/quark_l_h385_gT_dt_24_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5 ../F1M/corr/3pt_sp2/quark_l_h385_gT_dt_28_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5 ../F1M/corr/3pt_sp2/quark_l_h385_gT_dt_32_p2_0_ps2_0_g5P_g5P.fold.1835672416.h5'
OriginalTI='2 2 2 2 2'
OriginalTF='14 18 22 26 30'

##################################################

NumFields=words(FieldNames)

if( Save ne "" ) {
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  set output Save.'.pdf'
}

set key bottom right

fld=1
FieldName=word(FieldNames,fld)

# Get number of models
if( ModelMax < 0 ) {
  stats PlotFile using 'model' nooutput
  ModelMax=floor(STATS_max)
}
NumModels=ModelMax - ModelMin + 1
if( NumModels < 1 ) { print "Error: ".NumModels." models"; exit 1 }

# Work out ranges for each file
array FitTMin[NumModels]
array FitTMax[NumModels]
array FitSeqMin[NumModels]
array FitSeqMax[NumModels]
do for [model=1:NumModels] {
  stats PlotFile using (stringcolumn('field') eq FieldName && \
                        column('model') == (model-1+ModelMin) ? column('seq') : NaN ) \
                        : 't' nooutput
  FitSeqMin[model]=floor(STATS_min_x)
  FitSeqMax[model]=floor(STATS_max_x)
  FitTMin[model]=STATS_min_y
  FitTMax[model]=STATS_max_y
  #print "model=".(model-1+ModelMin).": seq=[".FitSeqMin[model].",".FitSeqMax[model]."], t=[".sprintf("%g",FitTMin[model]).",".sprintf("%g",FitTMax[model])."]"
}

set yrange [@my_yrange]
#set xrange [8.3:10.7]
#set xrange [7.8:11.2]
#set logscale y

# Work out layout
LMar=0.075
RMar=0.025
PlotWidth=1-LMar-RMar

# Total number of data points
NumDataPoints=0
do for [model=1:NumModels] { NumDataPoints=NumDataPoints+FitSeqMax[model] - FitSeqMin[model] + 1 }

array SubPlotWidth[NumModels]
array SubPlotLeft[NumModels]
array SubPlotRight[NumModels]
do for [model=1:NumModels] {
  SubPlotWidth[model]=PlotWidth * ( FitSeqMax[model] - FitSeqMin[model] + 1 ) / NumDataPoints
  SubPlotLeft[model]=model == 1 ? LMar : SubPlotRight[model - 1]
  SubPlotRight[model]=SubPlotLeft[model]+SubPlotWidth[model]
}

set multiplot layout 1, NumModels

do for [model=1:NumModels] {

if( MyTitle ne "" ) { set title word(MyTitle,model) }

set lmargin at screen SubPlotLeft[model]
set rmargin at screen SubPlotRight[model]

set xtics ceil(FitTMin[model]), 2, floor(FitTMax[model])

if( model == 2 ) { set format y ''; unset ylabel }

plot \
  PlotFile using (stringcolumn('field') eq word(FieldNames,fld) && \
                  column('model')==(model-1+ModelMin) ? column(xAxis) \
                  : NaN ) \
    :(column('theory_low')):(column('theory_high')) \
      with filledcurves notitle lc "skyblue" fs transparent solid 0.5\
  ,'' using (stringcolumn('field') eq word(FieldNames,fld) && \
              column('model')==(model-1+ModelMin) ? column(xAxis) : NaN ) \
          : (column('theory')) \
      with linespoints notitle lc "blue" dt 5 pt 4\
  ,'' using (stringcolumn('field') eq word(FieldNames,fld) && \
              column('model')==(model-1+ModelMin) ? column(xAxis) : NaN ) \
          : (column('data')):(column('data_low')):(column('data_high')) \
    with yerrorbars notitle lc "red" pt 13
}

unset multiplot
EOFMark
}

if [[ "$*" == "" ]];
then
  echo "$0"
  echo "Plot theory / data from fit."
  echo "Precede with optional modifiers (key=value):"
  echo "yrange    Vertical axis range (default: ALL DIFFERENT - i.e. BADLY LABELLED"
  echo "fields    Names of fields to display (default: log)"
  echo "fieldtext Names of fields fir legend (not used)"
  echo "title     Title for the plot"
  echo "size      of .pdf (default: 6in,3in)"
  echo "save      Filename to save (default: on screen)"
  echo "xAxis     Which field to show (default: t)"
  echo "mmin      Minimum model to show (default: 0)"
  echo "mmax      Maximum model to show (default: #last file)"
  exit 2
fi

PlotFile="$*"
PlotFunction
