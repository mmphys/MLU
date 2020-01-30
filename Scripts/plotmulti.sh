#!/bin/sh

# Loop through all the model files on the command-line performing plots
if [[ "$*" == "" ]]
then
  echo 'Specify a list of filenames to plot'
  echo 'optionally: field="im" to show imaginary part'
else

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark
#set term pdf

set xrange[${ti:=*}:${tf:=*}]
set yrange[${yrange:=*:*}]

FileNames="$*"
FieldName="${field:=y}"
MyUsing="\"t\":\"".FieldName."\":(column(\"".FieldName."\")-column(\"".FieldName."_low\")):(column(\"".FieldName."\")+column(\"".FieldName."_high\"))"
MySpec="using ".MyUsing." with yerrorbars"

RefVal=${ref:=-777}
if( RefVal != -777 ) {
  set arrow from graph 0, first RefVal to graph 1, first RefVal nohead front lc rgb "gray40" lw 0.25  dashtype "-"
  set label "E_0=".sprintf("%f",RefVal) at graph 0, first RefVal font "Arial,8" front textcolor "grey40" offset character 1.5, character 0.35
}

plot for [file in FileNames] file @MySpec title file noenhanced

EOFMark

fi
