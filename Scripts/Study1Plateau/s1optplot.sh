#!/opt/local/bin/bash

#Pick which "optimal" correlator we're plotting
PlotCorrName="${PlotCorr:=${mmplotdefault_corr}}_${ti:=${mmplotdefault_ti}}_${tf:=${mmplotdefault_tf}}"

#PlotDir="$HOME/data/201911HLFelix/Fit/One/"
#PlotBase="D"
#PlotCorr="uncorr"
#PlotSeed="2763836222"
#ti="4"
#tf="20"

times="${start:=0} ${step:=10} ${stop:=90}" #Must be evaluated outside of the subshell!!
times=`seq $times|tr '\n' ' '`
if (( !( $start <= 90 && $stop >= 90 || $start <= -90 && $stop >= -90 ) ))
then
  times="90 $times"
fi

###################################################
# Plot the masses of the 'optimal' correlator over a range
###################################################

gnuplot <<-EOFMark
# Begin Felix
#set term pdf size 5, 1.2
set key maxrows 3
#set yrange [0.99:1.035]
#set ytics 0.01
#set xrange [${ti}-3.2:${tf}+.2-2]
# End Felix

set term pdf #size 5, 1.2
set output "${mmplotdefault_base}.${PlotCorrName}.${start}_${stop}_${step}.theta.${mmplotdefault_seed}.pdf"
set ylabel '{/Helvetica:Italic aE}_{eff}'
set xlabel 't/a' offset 0,1

times="${times}"
numseries = words(times)
first_offset = ( numseries - 1 ) * 0.05 / 2

xr_min=${ti}+1-(first_offset+0.05); xr_max=${tf}-5+(first_offset+0.05); set xrange [xr_min:xr_max]
#set xrange [4.8:15.4]
set arrow from xr_min,0.99656 to xr_max,0.99656 nohead front lc rgb "gray40" lw 0.1 #lw 0.25  dashtype "-"

set pointsize 0.45
#set xlabel 'initial fit time'
#set palette defined ( 15 "blue", 16 "red", 17 "green")

#set title '"Optimal" operator [${PlotCorr}elated fit times ${ti} - ${tf}]'
#set ylabel 'Chi squared per degree of freedom'
#set logscale y
plot for [i=1:words(times)] \
 "${mmplotpath_model}${mmplotdefault_base}.${PlotCorrName}.theta_".word(times,i).".mass.$mmplotdefault_seed.$mmplotvar_dat" \
  using (\$1-first_offset+0.05*(i-1)):2:(\$2-\$3):(\$2+\$4) with yerrorbars title "θ=".word(times,i)
EOFMark
