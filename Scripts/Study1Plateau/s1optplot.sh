#!/bin/sh

# Default plot file
PlotDir="$HOME/data/201910Plat/fit/Z2/"
PlotBase="h1_l_p_0_0_0"
PlotCorr="uncorr"
PlotSeed="4147798751"
ti="6"
tf="16"

#echo "PlotDir=$PlotDir"
#echo "PlotBase=$PlotBase"
#echo "PlotCorr=$PlotCorr"
#echo "PlotSeed=$PlotSeed"
#echo "ti=$ti"
#echo "tf=$tf"

###################################################
# Plot the extracted energy levels and chi squared / d.o.f.
###################################################

gnuplot <<-EOFMark
#set term pdf

times="0 -90 -85 -80 -75"

set pointsize 0.6
#set xlabel 'initial fit time'
set xrange [${ti}-3.2:${tf}+.2-2]
#set palette defined ( 15 "blue", 16 "red", 17 "green")

#set output "${PlotBase}.${PlotCorr}.chisq.${PlotSeq}.pdf"
#ßset title '{/Times:Italic χ}^2 per d.o.f. dependence on initial/final ${PlotCorr}elated fit times ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
#set ylabel 'Chi squared per degree of freedom'
#set logscale y
#plot for [i=1:words(times)] "${PlotDir}${PlotBase}.${PlotCorr}_${ti}_${tf}.theta_".word(times,i).".mass.$PlotSeed.txt" \
  using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars title "theta=".word(times,i)
plot for [theta=-79:-73:1] "${PlotDir}${PlotBase}.${PlotCorr}_${ti}_${tf}.theta_".theta.".mass.$PlotSeed.txt" \
  using (\$1-0.15+0.05*(79+theta)):2:(\$2-\$3):(\$2+\$4) with yerrorbars title "theta=".theta
EOFMark
