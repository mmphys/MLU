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
local xAxis MyPlotField
for xAxis in ${x:-qSq EL}; do
for MyPlotField in ${PlotField:-data adjusted}; do
gnuplot <<-EOFMark

#Command-line options
PlotFile="$PlotFile"
my_xrange="${xrange:-*:*}"
my_yrange="${yrange:-*:*}"
my_ylabel="$YLabel"
yformat="$yformat"
FieldName="${field:-log}"
FieldNameText="${fieldtext:-effective mass}"
MyTitle="${title}"
PDFSize="${size:-6in,3in}"
Save="${save}"
xAxis="$xAxis"
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

if( Save ne "" ) {
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  OutFileName=Save.'_'.xAxis
  if( PlotField ne 'data' ) { OutFileName=OutFileName.'_'.PlotField }
  set output OutFileName.'.pdf'
}

if( RefText ne '' ) { RefText=RefText.', ' }
RefText=RefText.my_ylabel."(0)=$QSqZero, ".my_ylabel."(q^2_{max})=$QSqMax
if( PlotField ne 'data' ) { RefText=RefText.' ('.PlotField.')' }
set label 2 RefText at screen 0, screen 0 font ",10" front textcolor "blue" \
  offset character 0.5, 0.5

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

set xrange[@my_xrange]
set yrange[@my_yrange]
set pointintervalbox 0 # disables the white space around points in gnuplot 5.4.6 onwards
if( xAxis eq "EL" ) { set key top right } else { set key bottom right }
set xlabel xAxisName
set ylabel my_ylabel

# Styles I first used at Lattice 2022
#Ensembles="C1 C2 F1M M1 M2 M3"
Ensembles="C1 C2 F1M M1 M2 M3"
set linetype 1 ps 0.66 pt 6 lc rgb 0xFA60E4 #pink
set linetype 2 ps 0.66 pt 4 lc rgb 'dark-violet' # 0x9400D3
set linetype 3 ps 0.66 pt 14 lc rgb 0xE36C09 #orange
set linetype 4 ps 0.66 pt 10 lc rgb 0x00B050 #green
set linetype 5 ps 0.66 pt 12 lc rgb 0x0070C0 #blue
set linetype 6 ps 0.66 pt 8 lc rgb 0x003860 #dark navy
set linetype 7 ps 0.66 pt 7 lc rgb 0xFFC000 #yellow
set linetype 8 ps 0.66 pt 7 lc rgb 0xC00000 #red

#print "word(Ensembles,1)=".word(Ensembles,1)
#print "word(Ensembles,2)=".word(Ensembles,2)
#print "word(Ensembles,3)=".word(Ensembles,3)
plot PlotFile."_fit.txt" using (stringcolumn("field") eq xAxis ? column('x')*xScale : NaN) \
        :(column('y_low')):(column('y_high')) \
        with filledcurves notitle fc "skyblue" fs transparent solid 0.5, \
    for [i=1:6] PlotFile.".txt" \
        using (stringcolumn("ensemble") eq word(Ensembles,i) ? column(xAxis) * xScale : NaN) \
        :(column(PlotField)):(column(PlotField."_low")):(column(PlotField."_high")) \
        with yerrorbars title word(Ensembles,i) ls i, \
    PlotFile."_fit.txt" using (stringcolumn("field") eq xAxis ? column('x')*xScale : NaN) \
        :(column('y')) with lines title "Continuum" ls 8, \
    NaN with filledcurves title "Error" fc "skyblue" fs transparent solid 0.5
#    NaN with lines title 'Banana' lc rgb 0xC00000 bgnd "skyblue"

EOFMark
done
done
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

# If destination directory given: 1) Make it; and 2) ensure it ends in slash
if [[ -n $SaveDir ]]; then
  mkdir -p "$SaveDir"
  if [[ ${SaveDir: -1} != / ]]; then SaveDir="$SaveDir/"; fi
fi

# Save in specified directory, otherwise in same directory as source
# Input:  $1 filename to plot
# Output: save    Path to save to, without seed or extension (ie ending in .model)
function GetSaveName()
{
  local -n InPath=$1
  local Filename="${InPath##*/}"
  local IFS=.
  local -a NameParts
  read -ra NameParts <<< "$Filename"
  if (( ${#NameParts[@]} > 2 )); then
    NameParts[-2]=${NameParts[-2]%%_*}
    save="${NameParts[*]:0:${#NameParts[@]}-1}"
  else
    save="${NameParts[*]:0:0}.model"
  fi
  if [[ $Filename = $InPath ]]; then
    InPath="$save"
  else
    InPath="${InPath%/*}/$save"
  fi
  if [[ -v SaveDir || $Filename = $InPath ]]; then
    save="${SaveDir}$save"
  else
    save="${InPath%/*}/$save"
  fi
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

for PlotFile in "$@"
do
  if [[ -z "$PlotFile" ]]; then
    echo "Skipping empty filename"
  else
    if [[ ! -e "$PlotFile" && -e "$PlotFile.txt" ]]; then PlotFile+=.txt; fi
    if [ -v save ]; then PlotFile="${PlotFile%.*}"; else GetSaveName PlotFile; fi
    if [[ ! -e "$PlotFile.txt" ]]; then
      echo "Doesn't exist $PlotFile.txt"
    else
  (
    ff=${PlotFile##*/}
    ff=${ff#*.}
    ff=${ff%%.*}
    ff=${ff#*_}
    case $ff in
      fplus) YLabel="f_+";;
      f0) YLabel="f_0";;
    esac
    GetFitData
    PlotFunction
  ) fi
  fi
  unset save
done
