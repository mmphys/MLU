#!/usr/bin/env bash
. common_utility.sh

###################################################
# Plot theory - data from fit
#
# TODO: Fits with a single point don't show fit error range (filledcurves)
#
###################################################

PlotFunction()
{
gnuplot <<-EOFMark

#Command-line options
PlotFile="$PlotFile"
my_yrange="${yrange:-*:*}"
FieldName="${field:-log}"
FieldNameText="${fieldtext:-effective mass}"
MyTitle="${title}"
PDFSize="${size:-6in,3in}"
Save="${save}"
xAxis="${x:-t}"
ModelMin=${mmin:-0}
ModelMax=${mmax:--1}
OriginalTI="${ti}"
OriginalTF="${tf}"
OriginalFiles="${files}"
DoLog=0${log+1}

xOffset=FieldName eq "log" ? -0.5 : 0.0
NumTI=words( OriginalTI )
NumTF=words( OriginalTF )
#print 'NumTI='.NumTI.', OriginalTI='.OriginalTI
#print 'NumTF='.NumTF.', OriginalTF='.OriginalTF
#print OriginalFiles

if( Save ne "" ) {
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  set output Save.'.pdf'
}

set key bottom right

# Get number of models
if( ModelMax < 0 ) {
  stats PlotFile using 'model' nooutput
  ModelMax=floor(STATS_max)
}
NumModels=ModelMax - ModelMin + 1
if( NumModels < 1 ) { print "Error: ".NumModels." models"; exit 1 }

# This function sets the X-Axis to the requested field
GetUsingX(col)='using (stringcolumn("field") eq FieldName && column("model")==(model-1+ModelMin) ? column('.col.') : NaN )'

# Work out ranges for each file
array FitTMin[NumModels]
array FitTMax[NumModels]
array FitSeqMin[NumModels]
array FitSeqMax[NumModels]
array PlotTMin[NumModels]
array PlotTMax[NumModels]
array PlotTWidth[NumModels]
TotalPlotTWidth=0
UsingX=GetUsingX("'seq'")
do for [model=1:NumModels] {
  stats PlotFile @UsingX : 't' nooutput
  FitSeqMin[model]=floor(STATS_min_x)
  FitSeqMax[model]=floor(STATS_max_x)
  FitTMin[model]=STATS_min_y
  FitTMax[model]=STATS_max_y
  #print "model=".(model-1+ModelMin).": seq=[".FitSeqMin[model].",".FitSeqMax[model]."], t=[".sprintf("%g",FitTMin[model]).",".sprintf("%g",FitTMax[model])."]"
  PlotTMin[model]=(model <= NumTI) ? ( word(OriginalTI,model)+0 ) : ( FitTMin[model] - 0.5 )
  PlotTMax[model]=(model <= NumTF) ? ( word(OriginalTF,model)+0 ) : ( FitTMax[model] + 0.5 )
  #print 'PlotTMin[model]='.sprintf("%g",PlotTMin[model])
  #print 'PlotTMax[model]='.sprintf("%g",PlotTMax[model])
  PlotTWidth[model]=PlotTMax[model]-PlotTMin[model]
  #print 'PlotTWidth[model]='.sprintf("%g",PlotTWidth[model])
  TotalPlotTWidth=TotalPlotTWidth+PlotTWidth[model]
}
set yrange [@my_yrange]
#set xrange [8.3:10.7]
#set xrange [7.8:11.2]
if( DoLog ) { set logscale y }

# Work out layout
LMar=0.075
RMar=0.025
PlotWidth=1-LMar-RMar

array SubPlotWidth[NumModels]
array SubPlotLeft[NumModels]
array SubPlotRight[NumModels]
do for [model=1:NumModels] {
  SubPlotWidth[model]=PlotWidth * PlotTWidth[model] / TotalPlotTWidth
  SubPlotLeft[model]=model == 1 ? LMar : SubPlotRight[model - 1]
  SubPlotRight[model]=SubPlotLeft[model]+SubPlotWidth[model]
}

PointSizeData="0.6"
PointSizeTheory="0.72"
UsingX=GetUsingX("xAxis")

PlotPrefix='plot PlotFile '.UsingX
PlotPrefix=PlotPrefix.':(column("theory_low")):(column("theory_high"))'
PlotPrefix=PlotPrefix.' with filledcurves notitle lc "skyblue" fs transparent solid 0.5'
PlotPrefix=PlotPrefix.', "" '.UsingX
PlotPrefix=PlotPrefix.': (column("theory"))'
PlotPrefix=PlotPrefix.'with linespoints notitle lc "blue" dt 5 pt 4 ps '.PointSizeTheory

PlotSuffix=', PlotFile '.UsingX
PlotSuffix=PlotSuffix.': (column("data")):(column("data_low")):(column("data_high"))'
PlotSuffix=PlotSuffix.' with yerrorbars notitle lc "red" pt 13 ps '.PointSizeData

set multiplot layout 1, NumModels

do for [model=1:NumModels] {

if( MyTitle ne "" ) { set title word(MyTitle,model) }

set lmargin at screen SubPlotLeft[model]
set rmargin at screen SubPlotRight[model]

set xtics ceil(PlotTMin[model]), 2, floor(PlotTMax[model])

if( model == 2 ) { set format y ''; unset ylabel }

ThisFile=word( OriginalFiles, model )
if( ThisFile eq '' ) {
  ExtraPlot=''
} else {
  OtherFieldName=(FieldName eq "log") ? "exp" : FieldName
  ColumnX="(column('t')+(".sprintf("%g",xOffset)."))"
  #print 'ColumnX='.ColumnX
  ConditionA=ColumnX.'>='.sprintf("%g",PlotTMin[model]).' && '.ColumnX.'<'.sprintf("%g",FitTMin[model])
  #print 'ConditionA='.ConditionA
  ConditionB=ColumnX.'<='.sprintf("%g",PlotTMax[model]).' && '.ColumnX.'>'.sprintf("%g",FitTMax[model])
  #print 'ConditionB='.ConditionB
  Condition='(('.ConditionA.') || ('.ConditionB.'))'
  #print 'Condition='.Condition
  ExtraPlot=', "'.ThisFile.'" using ('.Condition.' ? '.ColumnX.' : NaN )'
  ExtraPlot=ExtraPlot.": (column('".OtherFieldName."'))"
  ExtraPlot=ExtraPlot.": (column('".OtherFieldName."_low'))"
  ExtraPlot=ExtraPlot.": (column('".OtherFieldName."_high'))"
  ExtraPlot=ExtraPlot." with yerrorbars notitle lc 'black' pt 13 ps ".PointSizeData
  #print 'ExtraPlot='.ExtraPlot
}

eval PlotPrefix.ExtraPlot.PlotSuffix
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
