#!/usr/bin/env bash
. common_utility.sh

#Pick which "optimal" correlator we're plotting
if [ "${fix+x}" == "x" ]; then FixOtherAngle=1; else FixOtherAngle=0; fi
DefaultAngle="${angle:-theta}"

# Loop through all the model files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ ( "$mmplotfile_type" == "bootstrap" || "$mmplotfile_type" == "fold" ) \
   && ( ${#mmplotfile_ops[@]} == 4 || ( ${#mmplotfile_ops[@]} == 6 && ( "${mmplotfile_ops[0]}" == "$DefaultAngle" \
                                                                     || "${mmplotfile_ops[2]}" == "$DefaultAngle" ) ) ) \
   && "$mmplotfile_ext" == "txt" ]]; then #Silently skip pdfs

# Make prefix and suffix by removing the angle we are varying
PlotTitle="${mmplotfile_base}.${mmplotfile_corr_all}"
PlotPrefix="${PlotTitle}."
PlotSuffix="_${mmplotfile_ops[-2]}_${mmplotfile_ops[-1]}"
unset OtherAngle
if [[ ${#mmplotfile_ops[@]} == 4 ]]; then
  angle="${mmplotfile_ops[0]}"
  AnglePrefix="${angle}"
else
  angle="$DefaultAngle"
  if [[ "$FixOtherAngle" == "0" ]]; then
    PlotTitle="${PlotTitle} ${mmplotfile_ops[0]}=${mmplotfile_ops[2]}"
    AnglePrefix="${mmplotfile_ops[0]}"
    OtherAngle="${mmplotfile_ops[2]}"
  else
    AnglePrefix="${angle}"
    if [[ "${mmplotfile_ops[0]}" == "$angle" ]]; then
      PlotSuffix="_${mmplotfile_ops[2]}_${mmplotfile_ops[3]}${PlotSuffix}"
      PlotTitle="${PlotTitle} ${mmplotfile_ops[2]}=${mmplotfile_ops[3]}"
    else
      PlotPrefix="${PlotPrefix}${mmplotfile_ops[0]}_${mmplotfile_ops[1]}_"
      PlotTitle="${PlotTitle} ${mmplotfile_ops[0]}=${mmplotfile_ops[1]}"
    fi
  fi
fi
PlotPrefix="${PlotPrefix}${AnglePrefix}_"
PlotSuffix="${PlotSuffix}."

if [[ "$LastFileCombo" == "${PlotPrefix}${PlotSuffix}${mmplotfile_type}${mmplotfile_seed}" ]]; then
  echo "Ignoring repetition of $PlotFile"
else
LastFileCombo="${PlotPrefix}${PlotSuffix}${mmplotfile_type}${mmplotfile_seed}"

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
OtherAngle="$OtherAngle"
angle="${angle}"
PlotTitle="$PlotTitle"
PlotPrefix="$PlotPrefix"
PlotSuffix="$PlotSuffix"
PlotType="${mmplotfile_type}"
PlotSeed=".${mmplotfile_seed}."
my_xrange="${ti:-*}:${tf:-*}"
my_yrange="${yrange:-*:*}"
my_xtics="${xtics}"
#my_key="${key:-bottom right}"
my_key="${key:-top right}"
times="${times}"
do_log=${log:-0}

angle_short="θ"
if( angle eq "phi" ) { angle_short="φ" }

# Begin Felix
#set term pdf size 5, 1.2
set key @my_key maxrows 3
#set yrange [0.99:1.035]
#set ytics 0.01
#set xrange [${ti}-3.2:${tf}+.2-2]
# End Felix

set term pdf #size 5, 1.2
FieldNameFile=FieldName
if( do_log ) {
set logscale y
FieldNameFile=FieldNameFile."_log"
}
OutFile=PlotPrefix
if( OtherAngle ne "" ) { OutFile=OutFile.OtherAngle."_" }
OutFile=OutFile."${start}_${stop}_${step}".PlotSuffix.FieldNameFile.PlotSeed."pdf"
set output OutFile
set ylabel '{/Helvetica:Italic aE}_{eff}'
set xlabel 't/a' offset 0,1
set yrange [@my_yrange]
if( my_xtics ne "" ) { set xtics @my_xtics }

numseries = words(times)
first_offset = ( numseries - 1 ) * 0.05 / 2
if( OtherAngle eq "" ) { TimeString=times } else {
  TimeString=""
  do for [i=1:words(times)] {
    ThisString=word(times,i)."_".OtherAngle."_".word(times,i)
    if( i == 1 ) { TimeString=ThisString } else { TimeString=TimeString." ".ThisString }
  }
}

set xrange [@my_xrange]
set object 1 rect from graph 0, first 0.99561 to graph 1, first 0.99751 fs solid 0.05 noborder fc rgb "gray10" behind
set arrow from graph 0, first 0.99656 to graph 1, first 0.99656 nohead front lc rgb "gray40" lw 0.25  dashtype "-"
#set label 2 "E_0=0.99656(95), ArXiv:1812.08791" at screen 0, 0 font "Arial,8" front textcolor "grey40" offset character 1, character 0.75

set pointsize 0.45
#set xlabel 'initial fit time'
#set palette defined ( 15 "blue", 16 "red", 17 "green")

set title FieldTitle." ".PlotTitle noenhanced
#set ylabel 'Chi squared per degree of freedom'
#set logscale y
FieldLow=FieldName."_low"
FieldHigh=FieldName."_high"
plot for [i=1:words(times)] \
 "${mmplotfile_path}".PlotPrefix.word(TimeString,i).PlotSuffix.PlotType.PlotSeed."$mmplotvar_txt" \
  using (column("t")-first_offset+0.05*(i-1)):(column(FieldName)):(column(FieldLow)):(column(FieldHigh)) with yerrorbars title angle_short."=".word(times,i)
EOFMark
fi;fi;done
