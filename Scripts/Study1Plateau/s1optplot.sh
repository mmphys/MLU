#!/bin/sh

# Default plot file
PlotDir="$HOME/data/201910Plat/fit/Z2/"
PlotBase="h1_l_p_0_0_0"
PlotCorr="uncorr"
PlotSeed="4147798751"
ti="6"
tf="16"

#PlotDir="$HOME/data/201911HLFelix/Fit/One/"
#PlotBase="D"
#PlotCorr="uncorr"
#PlotSeed="2763836222"
#ti="4"
#tf="20"

#echo "PlotDir=$PlotDir"
#echo "PlotBase=$PlotBase"
#echo "PlotCorr=$PlotCorr"
#echo "PlotSeed=$PlotSeed"
#echo "ti=$ti"
#echo "tf=$tf"

###################################################
# Plot the masses of the 'optimal' correlator over a range
###################################################

gnuplot <<-EOFMark
theta_start=90
theta_stop=130
theta_step=4
# For Felix
#set term pdf size 5, 1.2
#set key maxrows 4
#set yrange [0.995:1.02]
#set xrange [${ti}-3.2:${tf}+.2-2]
#For me
set term pdf
set xrange [${ti}-.2:${tf}+.2-2]
#set xrange [1:30]

set output "${PlotBase}.${PlotCorr}.theta_".theta_start."_".theta_stop."_".theta_step.".${PlotSeed}.pdf"

#times="0 -90 -85 -80 -75"

set pointsize 0.45
#set xlabel 'initial fit time'
#set palette defined ( 15 "blue", 16 "red", 17 "green")

#ßset title '{/Times:Italic χ}^2 per d.o.f. dependence on initial/final ${PlotCorr}elated fit times ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
#set ylabel 'Chi squared per degree of freedom'
#set logscale y
#plot for [i=1:words(times)] "${PlotDir}${PlotBase}.${PlotCorr}_${ti}_${tf}.theta_".word(times,i).".mass.$PlotSeed.txt" \
  using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars title "theta=".word(times,i)
plot for [theta=theta_start:theta_stop:theta_step] "${PlotDir}${PlotBase}.${PlotCorr}_${ti}_${tf}.theta_".theta.".mass.$PlotSeed.txt" \
  using (\$1-0.15+0.05*((theta-theta_start)/theta_step)):2:(\$2-\$3):(\$2+\$4) with yerrorbars title "theta=".theta
EOFMark
