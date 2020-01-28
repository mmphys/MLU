#!/bin/sh

if [[ "$*" != "" ]]
then
  echo 'Parameters unexpected'
else

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark
#set term pdf

set xrange[${xrange:=-0.2:4.2}]
set yrange[${yrange:=*:*}]
set key left

E0=0.99656
disp(x)=E0*E0+(2*pi/24)*(2*pi/24)*x

filename(n)="h1_l_p_".word("0_0_0 1_0_0 1_1_0 1_1_1 2_0_0",n+1).".corr.g5_gT5.params.4147798751.txt"

plot for [i=0:4] filename(i) index '[tf=19]' every ::4::4 using (i):(\$3*\$3):(column(-2)) linecolor variable title "p^2=".i, \
    disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion ?"


EOFMark

fi
