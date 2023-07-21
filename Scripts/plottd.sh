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
yformat="$yformat"
FieldName="${field:-log}"
FieldNameText="${fieldtext:-effective mass}"
MyTitle="${title}"
PDFSize="${size:-6in,3in}"
Save="${save}"
NameAdorn=${NameAdorn:-1}
xAxis="${x:-t}"
ModelMin=${mmin:-0}
ModelMax=${mmax:--1}
tExtra=${tExtra:-2}
OriginalTI="${ti}"
OriginalTF="${tf}"
ExtraFiles="${extra}"
DoLog=0${log+1}
RefVal="$RefVal"
RefText="$RefText"
Enhanced=1-0${NoEnh+1}

xOffset=FieldName eq "log" ? -0.5 : 0.0
NumTI=words( OriginalTI )
NumTF=words( OriginalTF )
#print 'NumTI='.NumTI.', OriginalTI='.OriginalTI
#print 'NumTF='.NumTF.', OriginalTF='.OriginalTF

# Get number of models
stats PlotFile using 'model' nooutput
MaxFileModel=floor(STATS_max)
if( ModelMax < 0 ) { ModelMax=MaxFileModel }
NumModels=ModelMax - ModelMin + 1
if( NumModels < 1 ) { print "Error: ".NumModels." models"; exit 1 }

if( Save ne "" ) {
  if( NameAdorn ) {
    Save=Save."_".FieldName
    if( ModelMin || ModelMax!=MaxFileModel ) { Save=Save."_".ModelMin."_".ModelMax }
    Save=Save.".$MLUSeed"
  }
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  set output Save.".pdf"
}

set pointintervalbox 0 # disables the white space around points in gnuplot 5.4.6 onwards
set key bottom right

# Returns Value if the condition is met, otherwise NaN
ValueIfCond(Value,cond)='(('.cond.') ? ('.Value.') : NaN)'
# Gets the value of a column if the condition is met, otherwise NaN
ColIfCond(col,cond)=ValueIfCond('column('.col.')',cond)
# Get a condition for whether the data point is in the fit
InFitPrefix='stringcolumn("field") eq FieldName && column("model")==(model-1+ModelMin)'
InFitCond(infit)=InFitPrefix.' && column("fitseq")'.(infit ? "!" : "=").'=-1'
JoinCond(A,Join,B)=(A eq '' || B eq '' ) ? A.B : '(('.A.') '.Join.' ('.B.'))'
# The second field is 1 or 0 to select points inside or outside the fit range
#GetSuffix(col,infit)=' && column("fitseq")'.(infit ? "!" : "=").'=-1 ? column('.col.') : NaN )'
#GetUsing(col,infit,cond)=GetPrefix.(cond eq "" ? "" : " && (".cond.")").GetSuffix(col,infit)
GetUsing(col,infit,cond)=ColIfCond(col,JoinCond(cond,"&&",InFitCond(infit)))

# Work out ranges for each file
array FitTMin[NumModels]
array FitTMax[NumModels]
array FitSeqMin[NumModels]
array FitSeqMax[NumModels]
array PlotTMin[NumModels]
array PlotTMax[NumModels]
array PlotTWidth[NumModels]
TotalPlotTWidth=0
UsingX='using '.GetUsing("'seq'",1,"")
do for [model=1:NumModels] {

  stats PlotFile @UsingX : 't' nooutput
  FitSeqMin[model]=floor(STATS_min_x)
  FitSeqMax[model]=floor(STATS_max_x)
  FitTMin[model]=STATS_min_y
  FitTMax[model]=STATS_max_y
  #print "model=".(model-1+ModelMin).": seq=[".FitSeqMin[model].",".FitSeqMax[model]."], t=[".sprintf("%g",FitTMin[model]).",".sprintf("%g",FitTMax[model])."]"
  PlotTMin[model]=((word(OriginalTI,model) ne "") ? (word(OriginalTI,model)+0) : \
                                  (FitTMin[model] > tExtra ? FitTMin[model] - tExtra : 0)) - 0.5
  PlotTMax[model]=((word(OriginalTF,model) ne "") ? (word(OriginalTF,model)+0) : \
                                  (FitTMax[model] + tExtra)) + 0.5
  #print 'PlotTMin[model]='.sprintf("%g",PlotTMin[model])
  #print 'PlotTMax[model]='.sprintf("%g",PlotTMax[model])
  PlotTWidth[model]=PlotTMax[model]-PlotTMin[model]
  #print 'PlotTWidth[model]='.sprintf("%g",PlotTWidth[model])
  TotalPlotTWidth=TotalPlotTWidth+PlotTWidth[model]
}
#print "TotalPlotTWidth=".sprintf("%g",TotalPlotTWidth)
set yrange [@my_yrange]
#set ytics font ', 8'
#set format y "%.6f"
if( yformat ne "" ) { set format y yformat }
#set xrange [8.3:10.7]
#set xrange [7.8:11.2]
if( DoLog ) { set logscale y }

# Work out layout
LMar=0.090
RMar=0.015
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
UsingX='using '.GetUsing("xAxis",1,"")
#print "UsingX=".UsingX

PlotPrefix='plot PlotFile '.UsingX
PlotPrefix=PlotPrefix.':(column("theory_low")):(column("theory_high"))'
PlotPrefix=PlotPrefix.' with filledcurves notitle lc "skyblue" fs transparent solid 0.5'
PlotPrefix=PlotPrefix.', "" '.UsingX
PlotPrefix=PlotPrefix.': (column("theory"))'
PlotPrefix=PlotPrefix.'with linespoints notitle lc "blue" dt 5 pt 4 ps '.PointSizeTheory

PlotSuffix=', PlotFile '.UsingX
PlotSuffix=PlotSuffix.':'.GetUsing('"data"',1,"")
PlotSuffix=PlotSuffix.':'.GetUsing('"data_low"',1,"")
PlotSuffix=PlotSuffix.':'.GetUsing('"data_high"',1,"")
PlotSuffix=PlotSuffix.' with yerrorbars notitle lc "red" pt 13 ps '.PointSizeData

if( RefText ne '' ) {
  set label 2 RefText at screen 0, screen 1 font ",10" front textcolor "blue" \
    offset character 0.5, -0.5
}

set multiplot layout 1, NumModels

if( RefVal ne '' ) {
  set object 1 rect from graph 0, first word(RefVal,4) to graph 1, first word(RefVal,6) \
      fs solid 1 noborder fc rgb 0xD0D0D0 behind
  set arrow from graph 0, first word(RefVal,5) to graph 1, first word(RefVal,5) \
      nohead front lc rgb "gray40" lw 0.25 dashtype "-"
}

#print "TotalPlotTWidth=".sprintf("%g",TotalPlotTWidth)

do for [model=1:NumModels] {

if( MyTitle ne "" ) {
  if( Enhanced ) {
    set title word(MyTitle,model)
  } else {
    set title word(MyTitle,model) noenhanced
  }
}

set lmargin at screen SubPlotLeft[model]
set rmargin at screen SubPlotRight[model]

set xrange [PlotTMin[model]:PlotTMax[model]]
set xtics ceil(PlotTMin[model]) + 1, 2, floor(PlotTMax[model])
if( TotalPlotTWidth >= 120 ) { set xtics font ", 7" }
else { if( TotalPlotTWidth >= 75 ) { set xtics font ", 9" } }

# Turn off y-axis and labels on second and subsequent models
if( model == 2 ) { set format y ''; unset ylabel }

OtherFieldName=(FieldName eq "log") ? "exp" : FieldName
ConditionTMin='>=PlotTMin[model]'
ConditionTMax='<=PlotTMax[model]'

# Plot data points outside fit range
Condition="column(xAxis)".ConditionTMin." && column(xAxis)".ConditionTMax
OriginalPlot=', PlotFile using '.GetUsing("xAxis",0,Condition)
OriginalPlot=OriginalPlot.':'.GetUsing('"data"',0,Condition)
OriginalPlot=OriginalPlot.':'.GetUsing('"data_low"',0,Condition)
OriginalPlot=OriginalPlot.':'.GetUsing('"data_high"',0,Condition)
OriginalPlot=OriginalPlot." with yerrorbars notitle lc 'black' pt 13 ps ".PointSizeData
#print 'OriginalPlot='.OriginalPlot

ThisFile=word( ExtraFiles, model )
if( ThisFile eq '' ) {
  ExtraPlot=''
} else {
  ColumnX="column('t')+xOffset"
  #print 'ColumnX='.ColumnX
  Condition='('.ColumnX.ConditionTMin.' && '.ColumnX.ConditionTMax.')'
  #print 'Condition='.Condition
  ExtraPlot=', "'.ThisFile.'" using '.ValueIfCond(ColumnX.'+0.1',Condition)
  ExtraPlot=ExtraPlot.": ".ColIfCond("OtherFieldName",Condition)
  ExtraPlot=ExtraPlot.": ".ColIfCond("OtherFieldName.'_low'",Condition)
  ExtraPlot=ExtraPlot.": ".ColIfCond("OtherFieldName.'_high'",Condition)
  ExtraPlot=ExtraPlot." with yerrorbars notitle lc 'gold' pt 6 ps ".PointSizeData
  #print 'ExtraPlot='.ExtraPlot
}

#print "model=".model
#print "PlotPrefix=".PlotPrefix
#print "ExtraPlot=".ExtraPlot
#print "OriginalPlot=".OriginalPlot
#print "PlotSuffix=".PlotSuffix
#print "========================================"

eval PlotPrefix.ExtraPlot.OriginalPlot.PlotSuffix
}

unset multiplot
EOFMark
}

unset bError
if [[ -v save && $# != 1 ]]; then
  bError=
  echo "Only 1 input file allowed when specifying save filename"
fi
if [[ -v save && -v SaveDir ]]; then
  bError=
  echo "Cannot specify both save and SaveDir"
fi

if [[ -v bError || -z $@ ]];
then
  echo "$0"
  echo "Plot theory / data from fit."
  echo "Precede with optional modifiers (key=value):"
  echo "yrange    Vertical axis range (default: ALL DIFFERENT - i.e. BADLY LABELLED"
  echo "yformat   Vertical axis format, e.g. %.5f for 5 decimal places"
  echo "log       Set to anything to make y-axis a log scale"
  echo "field     Name of field to display (default: log)"
  #echo "fieldtext Name of field for legend (not used)"
  echo "title     Title for each plot"
  echo "NoEnh     Set to anything to disable enhanced mode in titles"
  echo "size      of .pdf (default: 6in,3in)"
  echo "save      Filename to save (default: derived from PlotFile)"
  echo "SaveDir   Directory to save files to (default: cwd)"
  echo "xAxis     Which field to show (default: t)"
  echo "mmin      Minimum model to show (default: 0)"
  echo "mmax      Maximum model to show (default: #last file)"
  echo "tExtra    Number of extra data points to plot before and after fit (default 2)"
  echo "ti        List of t_min for each plot"
  echo "tf        List of t_max for each plot"
  echo "extra     List of extra files to plot - corresponding to each model"
  echo "RefVal    Reference value, as reported by GetColumn, see GetColumnValues()"
  echo "RefText   Reference text, see GetColumnValues()"
  exit 2
fi

# If destination directory given: 1) Make it; and 2) ensure it ends in slash
if [[ -n $SaveDir ]]; then
  mkdir -p "$SaveDir"
  if [[ ${SaveDir: -1} != / ]]; then SaveDir="$SaveDir/"; fi
fi

# Save in specified directory, otherwise in same directory as source
# Input:  $1 filename to plot
# Output: save    Path to save to, without seed or extension (ie ending in .model)
#         MLUSeed Seed from input filename (if present)
function GetSaveName()
{
  local InPath="$1"
  local Filename="${InPath##*/}"
  local IFS=.
  local -a NameParts
  read -ra NameParts <<< "$Filename"
  if (( ${#NameParts[@]} > 2 )); then MLUSeed="${NameParts[-2]}"; fi
  if (( ${#NameParts[@]} > 3 )); then
    if [[ ${NameParts[-3]: -3} = "_td" ]]; then NameParts[-3]=${NameParts[-3]:0:-3}; fi
    save="${NameParts[*]:0:${#NameParts[@]}-2}"
  else
    save="${NameParts[*]:0:1}.model"
  fi
  if [[ -v SaveDir || $Filename = $InPath ]]; then
    save="${SaveDir}$save"
  else
    save="${InPath%/*}/$save"
  fi
}

[ -v save ] && NameAdorn=0 || NameAdorn=1
for PlotFile in "$@"
do
  if [[ -z $PlotFile || ! -a $PlotFile ]]; then
    echo "Doesn't exist $PlotFile"
  else (
    if ! [ -v save ]; then GetSaveName "$PlotFile"; fi
    PlotFunction
  ) fi
done
