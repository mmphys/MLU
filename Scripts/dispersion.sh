#!/bin/sh

if [[ "$*" != "" ]]
then
  echo 'Parameters unexpected'
else

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark

tf="${tf:=19}"
q1="${q1:=h1}"
q2="${q2:=l}"
chi=${chi:=1.5}

E0=0.99656
L=24
kFactor=2 * pi / L
kFactorSq=kFactor * kFactor
#Momentum from k_x^2 + k_y^2 + k_z^2
#Momentum(KDotK)=kFactor * sqrt(KDotK)
#Dispersion relation from p
disp(KDotK)=E0 * E0 + kFactorSq * KDotK

filename(n)=q1."_".q2."_p_".word("0_0_0 1_0_0 1_1_0 1_1_1 2_0_0",n+1).".corr.g5_gT5.params.4147798751.txt"

set term pdf
set output "Disp_".tf.".pdf"

set xrange[${xrange:=-0.2:4.2}]
set yrange[${yrange:=*:*}]
set title q1." - ".q2." (a E_{eff})^2 vs k^2 for: {/Times:Italic χ}^2 ≤ ".sprintf("%g",chi)."; final fit time ".tf."; a E_0=".sprintf("%g",E0)
set ylabel '(a E_{eff})^2'
set xlabel 'k^2'
set xtics 1
set key top left
set pointsize 0.5

#plot for [i=0:4] filename(i) index "[tf=".tf."]" every ::4::4 using (i):(\$3*\$3):(column(-2)) linecolor variable title "k^2=".i, \
    disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2"

set label "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2" at graph 1, graph 0 font "Arial,12" front textcolor "grey40" offset character -1.5, character 1.5 right

plot for [i=0:4] filename(i) index "[tf=".tf."]" using (i):(column("ChiSqPerDof")<=chi?\$3*\$3:1/0):((\$3-\$4)*(\$3-\$4)):((\$3+\$5)*(\$3+\$5)):(\$2-9) with yerrorbars linecolor variable notitle, \
NaN ls 1 title "ti=10", \
NaN ls 2 title "ti=11", \
NaN ls 3 title "ti=12", \
NaN ls 4 title "ti=13", \
NaN ls 5 title "ti=14", \
NaN ls 6 title "ti=15", \
NaN ls 7 title "ti=16", \
NaN ls 8 title "ti=17", \
disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion"

EOFMark

fi
