#!/usr/bin/env bash
. common_utility.sh

#set -x
set -e

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
FileZ="$FileZ"
FileP="$FileP"
Save="$FileOut"
PDFSize="${size:-6in,3in}"
YField="${YField:-yNoPole}"
xAxis="${xAxis:-EL}"
my_xrange="${xrange:-*:*}"
my_yrange="${yrange:-*:*}"
my_ylabel="$YLabel"
yformat="$yformat"
MyTitle="${title}"
ModelMin=${mmin:-0}
ModelMax=${mmax:--1}
tExtra=${tExtra:-2}
OriginalTI="${ti}"
OriginalTF="${tf}"
ExtraFiles="${extra}"
RefVal="$RefVal"
RefText="$RefText"
PlotField="$MyPlotField"

Lambda=1e9
InvLambda=1e0 / Lambda; InvLambdaSq=InvLambda * InvLambda

if( xAxis eq "qSq" ) {
  xAxisName="q^2 / GeV^2"
  xScale=InvLambdaSq
} else {
  if( xAxis eq "EL" ) {
    xAxisName="E_{light} / GeV"
    xScale=InvLambda
  } else {
    xAxisName=xAxis
    xScale=1
  }
}

if( Save ne "" ) {
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  set output Save.'.pdf'
}

set pointintervalbox 0 # disables the white space around points in gnuplot 5.4.6 onwards
set key top left
set xlabel xAxisName
set ylabel "form factor x pole"

# Styles I first used at Lattice 2022
set linetype 1 ps 0.66 pt 6 lc rgb 0xFA60E4 #pink
set linetype 2 ps 0.66 pt 4 lc rgb 'dark-violet' # 0x9400D3
set linetype 3 ps 0.66 pt 14 lc rgb 0xE36C09 #orange
set linetype 4 ps 0.66 pt 10 lc rgb 0x00B050 #green
set linetype 5 ps 0.66 pt 12 lc rgb 0x0070C0 #blue
set linetype 6 ps 0.66 pt 8 lc rgb 0x003860 #dark navy
set linetype 7 ps 0.66 pt 7 lc rgb 0xFFC000 #yellow
set linetype 8 ps 0.66 pt 7 lc rgb 0xC00000 #red

plot FileZ using (stringcolumn("field") eq xAxis ? column('x')*xScale : NaN) \
        :(column(YField.'_low')):(column(YField.'_high')) \
        with filledcurves notitle fc "skyblue" fs transparent solid 0.5, \
    FileZ using (stringcolumn("field") eq xAxis ? column('x')*xScale : NaN):(column(YField)) \
        with lines title 'f_0' lc "skyblue", \
    FileP using (stringcolumn("field") eq xAxis ? column('x')*xScale : NaN) \
        :(column(YField.'_low')):(column(YField.'_high')) \
        with filledcurves notitle fc "red" fs transparent solid 0.5, \
    FileP using (stringcolumn("field") eq xAxis ? column('x')*xScale : NaN):(column(YField)) \
        with lines title 'f_+' lc "red"

EOFMark
}

# Save in specified directory, otherwise in same directory as source
# Input:  $1 filename to plot
# Output: save    Path to save to, without seed or extension (ie ending in .model)
function GetSaveName()
{
  local InPath=$1
  local SaveDir
  local InName=${InPath##*/}
  if [ "$InName" = "$1" ]; then
    SaveDir='.'
  else
    SaveDir="${InPath%/*}"
  fi
  FileZ=$SaveDir/F3_K_Ds.corr_f0.g5P_g5W.model_fit.txt
  FileP=$SaveDir/F3_K_Ds.corr_fplus.g5P_g5W.model_fit.txt
  FileOut=$SaveDir/F3_K_Ds.corr_f0_fplus.g5P_g5W.NoPole
  if [ -r "$FileZ" ] && [ -r "$FileP" ]; then return 0; fi
  echo "Files don't exist for $InPath"
  return 1
}

function GetFitData()
{
  local InPath=${PlotFile}.h5
  if ! [ -e $InPath ]; then
    InPath=${PlotFile/_f0/_f0_fplus}.h5
    if ! [ -e $InPath ]; then
      InPath=${PlotFile/_fplus/_f0_fplus}.h5
    fi
  fi
  if [ -e $InPath ]; then
    ColumnValues=($(GetColumn --exact ChiSqPerDof,pValue,pValueH $InPath))
    RefText="χ²/dof=${ColumnValues[0*8+4]} (p-H=${ColumnValues[2*8+4]}, p-χ²=${ColumnValues[1*8+4]})"
  else
    unset ColumnValues
    unset RefText
  fi
  local -a A
  A=($(gawk -e '/# qSq 0 / {print $4; next}; /# qSq / {print $4; exit}' ${PlotFile}_fit.txt))
  QSqZero="${A[0]}"
  QSqMax="${A[1]}"
  #echo $QSqZero
}

if [[ -z $@ ]];
then
  echo "$0"
  echo "Plot theory / data from fit."
  echo "Precede with optional modifiers (key=value):"
  echo "yrange    Vertical axis range (default: ALL DIFFERENT - i.e. BADLY LABELLED"
  echo "yformat   Vertical axis format, e.g. %.5f for 5 decimal places"
  echo "title     Title for the plot"
  echo "size      of .pdf (default: 6in,3in)"
  echo "save      Filename to save (default: derived from PlotFile)"
  echo "SaveDir   Directory to save files to (default: cwd)"
  echo "x         Which fields to show (default: qSq EL)"
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

for PlotFile in "$@"
do
  if GetSaveName "$PlotFile"; then
    PlotFunction
  fi
done
