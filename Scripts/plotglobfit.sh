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
my_xrange="${xrange:-*:*}"
my_yrange="${yrange:-*:*}"
yformat="$yformat"
FieldName="${field:-log}"
FieldNameText="${fieldtext:-effective mass}"
MyTitle="${title}"
PDFSize="${size:-6in,3in}"
Save="${save}"
xAxis="${x:-qSq}"
ModelMin=${mmin:-0}
ModelMax=${mmax:--1}
tExtra=${tExtra:-2}
OriginalTI="${ti}"
OriginalTF="${tf}"
ExtraFiles="${extra}"
RefVal="$RefVal"
RefText="$RefText"

if( Save ne "" ) {
  if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
  eval "set term pdfcairo dashed ".PDFSize
  set output Save.".pdf"
}

set xrange[@my_xrange]
set yrange[@my_yrange]
set pointintervalbox 0 # disables the white space around points in gnuplot 5.4.6 onwards
set key bottom right
set xlabel "q^2 / GeV^2"
set ylabel "$FormFactor"

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
plot for [i=1:3] PlotFile \
        using (stringcolumn("ensemble") eq word(Ensembles,i) ? column(xAxis) * 1e-18 : NaN) \
        :(column("data")):(column("data_low")):(column("data_high")) \
        with yerrorbars title word(Ensembles,i) ls i

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
  echo "title     Title for the plot"
  echo "size      of .pdf (default: 6in,3in)"
  echo "save      Filename to save (default: derived from PlotFile)"
  echo "SaveDir   Directory to save files to (default: cwd)"
  echo "xAxis     Which field to show (default: qSq)"
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

for PlotFile in "$@"
do
  if [[ -z $PlotFile || ! -a $PlotFile ]]; then
    echo "Doesn't exist $PlotFile"
  else (
    if ! [ -v save ]; then GetSaveName "$PlotFile"; fi
    FormFactor=${PlotFile##*/}
    FormFactor=${FormFactor#*.}
    FormFactor=${FormFactor%%.*}
    FormFactor=${FormFactor#*_}
    case $FormFactor in
      fplus) FormFactor="f_+";;
      f0) FormFactor="f_0";;
    esac
    PlotFunction
  ) fi
done
