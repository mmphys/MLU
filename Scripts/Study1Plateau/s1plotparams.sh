#!/usr/bin/env bash
. common_utility.sh

# Loop through all the params files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ "$mmplotfile_type" == "params_sort" ]]
then
  # Quietly skip the sorted version of parameter file
  mmplotfile_type="$mmplotfile_type"
elif [[ "$mmplotfile_type" != "params" ]]
then
  echo "Not params file: $mmplotfile_name"
else

if [[ "$stat" == "" ]]; then do_stat=0; else do_stat=1; fi
if [[ "$timin" == "" ]]; then do_timin=0; else do_timin=1; fi
if [[ "$timax" == "" ]]; then do_timax=0; else do_timax=1; fi
if [[ "$log" == "" ]]; then do_log=0; else let do_log=$log; fi
MyColumnHeadings=`awk '/^#/ {next}; {print \$0;exit}' $PlotFile`
MyColumnHeadingsNoUS="${MyColumnHeadings//_/\\\\_}"

###################################################
# Plot the extracted energy levels and test statistic associated with the model
###################################################

gnuplot <<-EOFMark

PlotFile="$PlotFile"
PlotFileSort="${mmplotfile_path}${mmplotfile_base}.${mmplotfile_type}_sort.${mmplotfile_seed}.${mmplotfile_ext}"
do_stat=${do_stat:-0}
stat_limit=${stat:-0}
do_timin=${do_timin}
TIMin=${timin:-0}
do_timax=${do_timax}
TIMax=${timax:-0}
my_key="${key:-bottom left maxrows 2}"
my_yrange="${yrange:-*:*}"
MyColumnHeadings="${MyColumnHeadings}"
MyColumnHeadingsNoUS="${MyColumnHeadingsNoUS}"
my_xtics="${xtics}"
do_log=${do_log}
xAxisName="${xaxis:-pvalueH pvalue ChiSqPerDof}"
RefVal="${ref}"
RefErr="${err}"
RefText="${reftext}"
RefSuffix="${refsuffix}"
do_label="${label:+x}"

# Work out how many fields there are per column by checking for absolute minimum
GotAbsMin=(system("awk '! /#/ {print (\$0 ~ /E0_min/) ? 1 : 0; exit}' ".word(PlotFile,1)) eq "0") ? 0 : 1
FieldsPerColumn=GotAbsMin ? 6 : 4
#print "FieldsPerColumn=".sprintf("%d",FieldsPerColumn)

# Find a column marked E0
MyColumnHeadings=system("awk '/^#/ {next}; {print \$0;exit}' ".PlotFile)
#print MyColumnHeadings
#print "MyField=".word(MyColumnHeadings,5)
FieldOffset=1
while( word(MyColumnHeadings,FieldOffset) ne "E0" ) {
  if( word(MyColumnHeadings,FieldOffset) eq "" ) { print "Can't find field E0"; exit gnuplot }
  FieldOffset=FieldOffset+1
}
#print "E0 field in column ".FieldOffset
#Work out how many exponentials there were in the fit
NumExp=1
#while( word(MyColumnHeadings,FieldOffset+NumExp*FieldsPerColumn) eq "E".NumExp ) { NumExp=NumExp+1 }
SearchField=1
while( word(MyColumnHeadings,FieldOffset+SearchField*FieldsPerColumn) ne "" ) {
  ThisField=word(MyColumnHeadings,FieldOffset+SearchField*FieldsPerColumn)
  if( ThisField[1:1] eq "E" ) {
    ThisNumExp=ThisField[2:*]+1
    if( NumExp < ThisNumExp ) { NumExp = ThisNumExp }
  }
  SearchField=SearchField+1
}
#print "NumExp=".NumExp

# Now look for x-axis fields. Prefer fields at the start
xAxisNum=words(xAxisName) + 1
xAxisOffset=1
xAxisCol=-1
while( xAxisNum > 1 && word(MyColumnHeadings,xAxisOffset) ne "" ) {
  i = 1
  while( word(xAxisName,i) ne "" && word(xAxisName,i) ne word(MyColumnHeadings,xAxisOffset) ) {i=i+1}
  if( word(xAxisName,i) ne "" ) { xAxisNum=i; xAxisCol=xAxisOffset }
  xAxisOffset=xAxisOffset+1
}
if( xAxisCol==-1 ) {
  print "Error: couldn't find statistic column in: ".xAxisName
  exit gnuplot
}
xAxisName=word(xAxisName,xAxisNum)
#print "xAxisName=".xAxisName.", xAxisNum=".xAxisNum.", xAxisCol=".xAxisCol
if( xAxisName eq "ChiSqPerDof" ) {
  sStatDescr='{/Times:Italic χ}^2 per d.o.f.'
  sStatCond='>'
  sStatCondHuman='≤'
} else {
  sStatDescr=xAxisName
  sStatCond='<'
  sStatCondHuman='≥'
  if( do_stat == 0 ) { stat_limit=0.05; do_stat=1 }
}

if( do_stat ) { sStatDescr=sStatDescr.' ('.sStatCondHuman.' '.sprintf("%g",stat_limit).")" }
#FitName=' from '.NumExp."-exponential ${mmplotfile_corr}elated fit using ${mmplotfile_ops_all//_/ }"
FitName=' from '.NumExp."-exponential ${mmplotfile_corr}elated fit of ${mmplotfile_base}"
OutBase="${mmplotfile_base}.${mmplotfile_corr_all}.${mmplotfile_ops_all}."
OutSuffix="_".xAxisName.".${mmplotfile_seed}.pdf"

set term pdf
set pointsize 0.6
set xlabel sStatDescr
#set xrange [*:${stat:-*}]
#set palette defined ( 15 "blue", 16 "red", 17 "green")
set key @my_key
if( do_log ) { set logscale x }
if( my_xtics ne "" ) { set xtics @my_xtics }

stats PlotFile using "ti" nooutput
NBlock=STATS_blocks

set yrange[@my_yrange]

Condition=''
AppendCondition( s1, s2 )=s1.( strlen(s1) ? ' || ' : '' ).s2
if( do_timin ) { Condition=AppendCondition(Condition, 'column("ti") < TIMin' ) }
if( do_timax ) { Condition=AppendCondition(Condition, 'column("ti") > TIMax' ) }
if( do_stat ) { Condition=AppendCondition(Condition, 'column(xAxisCol) '.sStatCond.' stat_limit' ) }
if( strlen( Condition ) ) { Condition=Condition.' ? NaN : ' }

WithLabels='with labels font "Arial,6" rotate noenhanced left offset char 0, 0.25'

# Loop through all fields, plotting them

if( do_label eq "x" ) {
  LabelStart=NBlock
} else {
  LabelStart=0
}

while( word(MyColumnHeadings,FieldOffset) ne "" ) {
  MyField=word(MyColumnHeadings,FieldOffset)
  MyFieldNoUS=word(MyColumnHeadingsNoUS,FieldOffset)
  OutFile=OutBase.MyField.OutSuffix
  #print "MyField=".MyField.", OutFile=".OutFile
  set output OutFile
  #set label 1 OutFile noenhanced at screen 1, 0 right font "Arial,8" front textcolor "grey40" offset character -1, character 0.75
  set label 1 OutFile noenhanced at screen 1, 0.5 center rotate by -90 font "Arial,8" front textcolor "grey40" offset character -1.5, character 0
  IsEnergy=0
  if (MyField[1:1] eq "E") {
    EnergyLevel=MyField[2:*] + 0
    if( EnergyLevel < words(RefVal) ) {
      IsEnergy=1
      ThisRef=word(RefVal,EnergyLevel+1)+0
      ThisErr=word(RefErr,EnergyLevel+1)+0
      if( EnergyLevel < words(RefText) ) { ThisRefText=word(RefText,EnergyLevel+1)
      } else { ThisRefText=" ± ".word(RefErr,EnergyLevel+1) }
      ThisRefText=word(RefVal,EnergyLevel+1).ThisRefText.RefSuffix
      set object 1 rect from graph 0, first ThisRef-ThisErr to graph 1, first ThisRef+ThisErr \
        fs solid 0.05 noborder fc rgb "gray10" behind
      set arrow from graph 0, first ThisRef to graph 1, first ThisRef nohead front lc rgb "gray40" lw 0.25  dashtype "-"
      set label 2 "a ".MyField." = ".ThisRefText \
        at screen 0, 0 font "Arial,8" front textcolor "grey40" offset character 1, character 0.75
    }
  }
  MyTitle = MyFieldNoUS.FitName
  set title MyTitle noenhanced
  set ylabel MyFieldNoUS
  
  plot \
    for [idx=0:NBlock-1] PlotFile index idx using \
      (@Condition column(xAxisCol)):(@Condition column(MyField)):(@Condition column(MyField."_low")):(@Condition column(MyField."_high")) \
      with yerrorbars title columnheader(1), \
    for [idx=LabelStart:NBlock-1] '' index idx using (@Condition column(xAxisCol)):(@Condition column(MyField."_high")):(column("tfLabel")) \
      @WithLabels notitle

  if (IsEnergy) { unset object 1; unset arrow; unset label 2 }
  FieldOffset=FieldOffset + FieldsPerColumn
}

unset xrange
unset logscale x
#set yrange [*:${stat:-*}]
unset yrange
OutFile=OutBase.xAxisName.'-seq'.OutSuffix
set output OutFile
set label 1 OutFile noenhanced at screen 1, 0 right font "Arial,8" front textcolor "grey40" offset character -1, character 0.75
set title sStatDescr." ".FitName noenhanced
set ylabel sStatDescr
set xlabel 'Sequence'
set key top left
unset xtics

MyField=xAxisName
plot \
  for [idx=0:NBlock-1] PlotFile index idx using \
    (@Condition column(1)):(@Condition column(MyField)):(@Condition column(MyField."_low")):(@Condition column(MyField."_high")) \
    with yerrorbars title columnheader(1), \
  for [idx=LabelStart:NBlock-1] '' index idx using \
    (@Condition column(1)):(@Condition column(MyField."_high")):(column("tfLabel")) \
    @WithLabels notitle
EOFMark

fi
done
