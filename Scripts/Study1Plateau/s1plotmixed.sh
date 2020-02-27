#!/opt/local/bin/bash

#Pick which "optimal" correlator we're plotting

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ "$mmplotfile_type" == "bootstrap" && "${mmplotfile_ops_all:0:6}" == "theta_" && "$mmplotfile_ext" == "txt" ]]; then #Silently skip pdfs
if [[ "$LastFileCombo" == "${mmplotfile_base}.${mmplotfile_corr_all}..${mmplotfile_seed}." ]]; then
  echo "Ignoring repetition of $PlotFile"
else
PlotPrefix="${mmplotfile_base}.${mmplotfile_corr_all}."
PlotSuffix=".${mmplotfile_seed}."
LastFileCombo="${PlotPrefix}${PlotSuffix}"

times="${start:=0} ${step:=10} ${stop:=90}" #Must be evaluated outside of the subshell!!
times=`seq $times|tr '\n' ' '`
if (( !( $start <= 90 && $stop >= 90 || $start <= -90 && $stop >= -90 ) ))
then
  times="90 0 $times"
fi

case ${field:=cosh} in
  cosh | exp) field_title="$field mass";;
           *) field_title="$field";;
esac

###################################################
# Plot the masses of the 'optimal' correlator over a range
###################################################

gnuplot <<-EOFMark

FieldName="$field"
FieldTitle="$field_title"
PlotPrefix="$PlotPrefix"
PlotSuffix="$PlotSuffix"
my_xrange="${ti:-*}:${tf:-*}"
#my_key="${key:-bottom right}"
my_key="${key:-top right}"
times="${times}"

# Begin Felix
#set term pdf size 5, 1.2
set key @my_key maxrows 3
#set yrange [0.99:1.035]
#set ytics 0.01
#set xrange [${ti}-3.2:${tf}+.2-2]
# End Felix

set term pdf #size 5, 1.2
set output PlotPrefix."theta.${start}_${stop}_${step}.".FieldName.PlotSuffix."pdf"
set ylabel '{/Helvetica:Italic aE}_{eff}'
set xlabel 't/a' offset 0,1

numseries = words(times)
first_offset = ( numseries - 1 ) * 0.05 / 2

set xrange [@my_xrange]
set arrow from graph 0, first 0.99656 to graph 1, first 0.99656 nohead front lc rgb "gray40" lw 0.25  dashtype "-"
set label 2 "E_0=0.99656(95), ArXiv:1812.08791" at screen 0, 0 font "Arial,8" front textcolor "grey40" offset character 1, character 0.75

set pointsize 0.45
#set xlabel 'initial fit time'
#set palette defined ( 15 "blue", 16 "red", 17 "green")

set title FieldTitle." from mixed operator ".PlotPrefix noenhanced
#set ylabel 'Chi squared per degree of freedom'
#set logscale y
FieldLow=FieldName."_low"
FieldHigh=FieldName."_high"
print "${mmplotfile_path}".PlotPrefix."theta_".word(times,1).".bootstrap".PlotSuffix."$mmplotvar_dat"
plot for [i=1:words(times)] \
 "${mmplotfile_path}".PlotPrefix."theta_".word(times,i).".bootstrap".PlotSuffix."$mmplotvar_dat" \
  using (column("t")-first_offset+0.05*(i-1)):(column(FieldName)):(column(FieldLow)):(column(FieldHigh)) with yerrorbars title "θ=".word(times,i)
EOFMark
fi;fi;done
