#!/usr/bin/env bash
. common_utility.sh

#Optional environment variables to set before running
# corr2: set some useful defaults for two correlators

# Important defaults here
min=${min:-5}
max=${max:-15}
tione=${tione:-$min}
key=${key:-bottom left maxrows 2}
series=${series:-ti}
size=${size:-"7in,2in"}
if [ "$serieslbl" == "" ]
then
  case "$series" in
    ti)  serieslbl="tf";;
    tf)  serieslbl="ti";;
    ti1) serieslbl="tf1";;
    tf1) serieslbl="ti1";;
    *)   serieslbl="tf";;
  esac
fi
# Defaults if we ask for two-correlator defaults
if [ -v corr2 ]; then
  [ -v where ] || where='column("ti") == column("ti1")'
  [ -v extralbl ] || extralbl='column("ti") == column("ti1") ? "" : "(".stringcolumn("ti1").")"'
fi

if [ ${#@} == 0 ]; then
  echo "$0"
  echo "Plot fit results (i.e. parameters from model)"
  echo "Precede with optional modifiers (key=value):"
  echo "sort      options to pass to sort utility to sort before plot"
  echo "          e.g. to sort on columns 2 through 5, specify '-V -k 2,5'"
  echo "min       Minimum TI ($min)"
  echo "max       Maximum TI ($max)"
  echo "tione     Which colour should be number 1 ($tione)"
  echo "series    Which field is the series category ($series)"
  echo "serieslbl Which field is the series label ($serieslbl)"
  echo "extralbl  Extra fields to plot in label ($extralbl)"
  echo "where     Only plot points where constraint met ($where)"
  echo "key       Defaults for legend ($key)"
  echo "seq       Add sequence number to label"
  echo "pvalue    Add pvalue to point label"
  echo "size      Size of pdf ($size)"
  echo "save      Set to save files. If not empty, this is suffix for filename($save)"
fi

# Loop through all the params files on the command-line performing plots
for PlotFile; do
PlotPathSplit "$PlotFile" 6
if [[ "${mmplotfile_type:0:6}" != "params" ]]
then
  echo "Expected 'params*' file not '$mmplotfile_type'"
else

MyColumnHeadings=`awk '/^#/ {next}; {print \$0;exit}' $PlotFile`
MyColumnHeadingsNoUS="${MyColumnHeadings//_/\\\\_}"

# sort the file if requested
if [ -v sort ]
then
  # Make a temporary file for the header
  TmpFile=$(mktemp tmp.Header.XXXXXX)
  awk '/^#/ {print;next};{if(c++){exit}print}' $PlotFile > $TmpFile
  # Now append the sorted contents
  awk '/^ti=/ {next}; /^#/ {next}; !z {z=1}; {print}' $PlotFile | sort $sort >> $TmpFile
  # Now select the temporary file
  PlotFile=$TmpFile
fi

###################################################
# Plot the extracted energy levels and test statistic associated with the model
###################################################

gnuplot -p <<-EOFMark

PlotFile="$PlotFile"
PlotFileSort="${mmplotfile_path}${mmplotfile_base}.${mmplotfile_type}_sort.${mmplotfile_seed}.${mmplotfile_ext}"
my_key="${key:-bottom left maxrows 2}"
MyColumnHeadings="${MyColumnHeadings}"
MyColumnHeadingsNoUS="${MyColumnHeadingsNoUS}"
my_xtics="$xtics"
PDFSize="$size"
MyInteractive=1-0${save+1}
MySeries="${series}"
MySeriesMin=${min}
MySeriesMax=${max}
MyFirstColour=${tione}

OutBase="${mmplotfile_base}.${mmplotfile_corr_all}.${mmplotfile_ops_all}.${save:+${save}_}"
#if( 0${save:+1} ) { OutSuffix="$save_".OutSuffix }
OutSuffix=".${mmplotfile_seed}.pdf"

# Find the column number of the first of the value fields
FieldNames="E0"
FieldOffset=1
while( word(MyColumnHeadings,FieldOffset) ne FieldNames ) {
  if( word(MyColumnHeadings,FieldOffset) eq "" ) { print "Can't find field ".FieldNames; exit gnuplot }
  FieldOffset=FieldOffset+1
}
#print FieldNames." is field ".FieldOffset

# Work out how many fields there are per column by checking for absolute minimum
FieldsPerColumn=word(MyColumnHeadings,FieldOffset-2) eq "E0_min" ? 6 : 4
#print "FieldsPerColumn=".FieldsPerColumn

# Get a list of all the available fields
NumFields=words(FieldNames)
while( word(MyColumnHeadings,NumFields*FieldsPerColumn+FieldOffset) ne "ChiSqPerDof" && \
       word(MyColumnHeadings,NumFields*FieldsPerColumn+FieldOffset) ne "" ) {
  FieldNames=FieldNames." ".word(MyColumnHeadings,FieldOffset+NumFields*FieldsPerColumn)
  NumFields=NumFields+1
}
#print "NumFields=".NumFields
#print "FieldNames=".FieldNames

Condition="column('".MySeries."') == idx"
if( 0${where:+1} ) { Condition=Condition.' && ($where)' }
#print "Condition=".Condition
ConditionLong='! ('.Condition.') ? NaN :'
WithLabels='with labels font ",5" rotate noenhanced left offset char 0, 0.25'
XPos='idx+(Count[MyIndex(idx)] <= 3 ? 0.4 : 0.8)/(Count[MyIndex(idx)]-1)*(Seq[MyIndex(idx)] - (Count[MyIndex(idx)]+1)*0.5)'
MyLabelText='stringcolumn("${serieslbl}")'
if( 0${extralbl:+1} ) { MyLabelText=MyLabelText.'.( ${extralbl} )' }
if( 0${pvalue+1} ) { MyLabelText=MyLabelText.'." p=".stringcolumn("pvalueH")' }
if( 0${seq+1} ) { MyLabelText=MyLabelText.'."  ".Seq[MyIndex(idx)]."/".Count[MyIndex(idx)]' }
#print "MyLabelText=".MyLabelText

# I'm going to need a sequence number and the total count within each series
MyIndex(idx) = idx - MySeriesMin + 1
MySeriesCount=MyIndex(MySeriesMax)
array Count[MySeriesCount]
array Seq[MySeriesCount]

set key $key
if( "$xrange" ne "" ) { set xrange [${xrange}] }
if( "$yrange" ne "" ) { set yrange [${yrange}] }

# Default line width and point type
DefWidth=0.75
DefType=7

# Default point size
if( MyInteractive ) { DefSize=0.85 } else { DefSize=0.4 }

# Looking at one file interactively
#MyFieldNum=1
#while( word(FieldNames,MyFieldNum) ne MyField ) {
  #if( word(MyColumnHeadings,MyFieldNum) eq "" ) { print "Can't find field ".MyFieldNum; exit gnuplot }
  #MyFieldNum=MyFieldNum+1
#}

# Looping over all files and saving them
do for [MyFieldNum=NumFields:1:-1] {
MyColumn=(MyFieldNum-1)*FieldsPerColumn+FieldOffset
MyField = word(FieldNames,MyFieldNum)
#print MyField." is field ".MyFieldNum.", column number ".MyColumn

if( MyInteractive ) {
  # Looking at one file interactively
  set terminal wxt MyFieldNum - 1 title MyField
  #set title MyField
} else {
  # Looping over all files and saving them
  SetTerm='set term pdfcairo font "Arial,12"'
  if( PDFSize ne "" ) { SetTerm=SetTerm."size ".PDFSize }
  #print SetTerm
  eval SetTerm
  OutFile=OutBase.MyField.OutSuffix
  set output OutFile
}

plot \
  for [idx=MySeriesMin:MySeriesMax] PlotFile using \
    (Count[MyIndex(idx)]=0, Seq[MyIndex(idx)]=0, NaN) : (NaN) notitle, \
  for [idx=MySeriesMin:MySeriesMax] PlotFile using \
    ((@Condition) ? Count[MyIndex(idx)]=Count[MyIndex(idx)] + 1 : 0, NaN) : (NaN) notitle, \
  for [idx=MySeriesMin:MySeriesMax] '' using \
    ((@Condition) ? Seq[MyIndex(idx)]=Seq[MyIndex(idx)] + 1 : 0, @ConditionLong @XPos) : \
    (column(MyColumn)) : (column(MyColumn-1)) : (column(MyColumn+1)) \
    with yerrorbars \
    linestyle idx - MyFirstColour + 1 \
    lw DefWidth pt DefType ps DefSize \
    title "${series}=".idx, \
  for [idx=MySeriesMin:MySeriesMax] '' using (Seq[MyIndex(idx)]=0, NaN) : (NaN) notitle, \
  for [idx=MySeriesMin:MySeriesMax] '' using \
    ((@Condition) ? Seq[MyIndex(idx)]=Seq[MyIndex(idx)] + 1 : 0, @ConditionLong @XPos) : \
    (column(MyColumn+1)) : (@MyLabelText) \
    @WithLabels notitle

}

if( MyInteractive ) { pause mouse close }
EOFMark

# Kill my temporary file
[ -v sort ] && rm $TmpFile

fi
done
