#!/usr/bin/env bash

#set -x
set -e

###################################################

# User input

###################################################

###################################################

# Plot correlation matrix

###################################################

function DoPlot()
{
local InFile="$1"
gnuplot <<-EOFMark

InFile="$InFile"
Save="${InFile%.txt}.pdf"
PDFSize="${size:-5in,4in}"

if( PDFSize ne "" ) { PDFSize="size ".PDFSize }
eval "set term pdfcairo ".PDFSize
set output Save

set xtics rotate noenhanced font 'Arial,8'
set ytics noenhanced font 'Arial,8'
#set title noenhanced "${InFile##*/}";
unset key
set cbrange[-1:1]
plot InFile matrix columnheaders rowheaders with image pixels

EOFMark
}

###################################################

# Main loop

###################################################

for PlotFile in $@
do
  DoPlot "$PlotFile"
done

