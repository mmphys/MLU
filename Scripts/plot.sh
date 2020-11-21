#!/usr/bin/env bash
. common_utility.sh

###################################################
# Make a plot of the forward and backward-propagating waves
###################################################

PlotFunction()
{
gnuplot <<-EOFMark

#Command-line options
nt=${nt:-0}
do_title="${do_title:-0}"
my_title="${title}"
my_x_axis="${x:-column(1)}"
xargs="${xargs}"
my_xrange="${ti:-*}:${tf:-*}"
my_yrange="${yrange:-*:*}"
my_cbrange="${cbrange:--1:1}"
my_ylabel="${ylabel}"
my_xlabel="${xlabel}"
my_key="${key:-top right}"
FieldNames="${fields:-cosh}"
RefVal=${ref:--777}
RefErr=${err:--777}
RefText="${reftext}"
do_log=${log:-0}
SaveFile=${SaveFile:-0}
SaveFileName="$SaveFileName"
do_offset=${offset:-0.05}
do_whisker=${whisker:-0}

# Which field names do we plot?
f = 1
FieldNamesShort=""
while (f <= words(FieldNames) ) {
  FieldNamesShort=FieldNamesShort."_".word(FieldNames,f)
  f=f+1
}

# Are we plotting forward and / or backward waves?
fb_min=0
fb_max=0
array fb_prefix[2]
fb_prefix[1]=""
FBShort=""
FBText=""
if ( nt != 0) {
  # We want the backward propagating wave
  fb_max = 1
  fb_prefix[2]="back "
  if ( nt < 0 ) {
    # We ONLY want the backward propagating wave
    nt = -nt
    fb_min = 1
    FBShort="_b"
    FBText=" backward"
  } else {
    fb_prefix[1]="fwd "
    FBShort="_fb"
    FBText=" forward/backward"
  }
}

# Are we doing a log plot?
LogShort=""
if( do_log ) {
  LogShort="_log"
}

PlotFile="$PlotFile"
NumFiles=words(PlotFile)
FileType="${mmplotfile_type}"

# Work out whether we have absolute minimum ... and can therefore do box and whiskers
if( do_whisker ) {
do_whisker=(system("awk '! /#/ {print (\$0 ~ /".word(FieldNames,1)."_min/) ? 1 : 0; exit}' ".word(PlotFile,1)) eq "0") ? 0 : 1 }

# Decide on a title and output filename
if( SaveFile == 2 ) {
  MyTitle = SaveFileName
  SaveFileName=SaveFileName.".pdf"
} else {
  # We are only processing one file
  MyTitle="${mmplotfile_name_no_ext}".FBText
  SaveFileName="${mmplotfile_base}.".FileType
  if( FileType ne "cormat" ) { SaveFileName=SaveFileName.FieldNamesShort.LogShort.FBShort }
  SaveFileName=SaveFileName.".${mmplotfile_seed}.pdf"
}

if( SaveFile ) {
  if( NumFiles==1 && FileType eq "cormat" ) {
    set term pdf size 7.25,7 fontscale 0.4
  } else {
    set term pdf
  }
  set pointsize 0.5
  set output SaveFileName
  set label 1 SaveFileName noenhanced at screen 1, 0.5 center rotate by -90 font "Arial,8" front textcolor "grey40" offset character -1.3, 0
}

set key font "Arial,8" @my_key noenhanced
set xrange[@my_xrange]
set yrange[@my_yrange]
if( do_title ) {
  # This title might be set to empty string, meaning - no title
  if( my_title ne "" ) { set title my_title }
} else {
  set title MyTitle noenhanced
}
if( my_xlabel ne "" ) { set xlabel my_xlabel }
if( my_ylabel ne "" ) { set ylabel my_ylabel }

if( RefVal != -777 ) {
  RefErrString=""
  if( RefErr != -777 ) {
    set object 1 rect from graph 0, first RefVal - RefErr to graph 1, first RefVal + RefErr fs solid 0.05 noborder fc rgb "gray10" behind
    RefErrString=sprintf(" (%g)", RefErr)
  }
  set arrow from graph 0, first RefVal to graph 1, first RefVal nohead front lc rgb "gray40" lw 0.25  dashtype "-"
  if( RefText eq "" ) { RefText="Ref: ".sprintf("%g", RefVal).RefErrString }
  set label 2 RefText at screen 0, screen 0 font "Arial,8" front textcolor "grey40" offset character 0.5, 0.25
}

#AbsMin(y,low,high)=sgn(y) < 0 ? -(high) : low
#AbsMax(y,low,high)=sgn(y) < 0 ? -(low) : high

if( SaveFile != 2 ) {
  # Don't do this for multi-plots, as probably best to stick to colour scheme
  set linetype 1 lc rgb 'blue'
  set linetype 2 lc rgb 'red'
}

# The logging variable becomes a log
if( do_log ) { set logscale y }

# Work out how much to offset each series by
NumFields=words(FieldNames)
NumFB=fb_max - fb_min + 1
XF1=do_offset
XF2=(1 - NumFields*NumFB*NumFiles)/2*XF1

PlotWith="yerrorbars"
PlotUsing=""
if( NumFiles==1 && FileType eq "cormat" ) {
  set xtics rotate noenhanced
  set ytics noenhanced
  set cbrange [@my_cbrange]
  plot PlotFile matrix columnheaders rowheaders with image pixels
} else {
  if( do_log ) {
    PlotUsing='(abs(column(word(FieldNames,fld)))) : (sgn(column(word(FieldNames,fld)))*column(word(FieldNames,fld).( sgn(column(word(FieldNames,fld))) ? "_low" : "_high" ))) : (sgn(column(word(FieldNames,fld)))*column(word(FieldNames,fld).( sgn(column(word(FieldNames,fld))) ? "_high" : "_low" )))'
  } else {
    if( !do_whisker ) {
      PlotUsing='(column(word(FieldNames,fld))) : (column(word(FieldNames,fld)."_low")) : (column(word(FieldNames,fld)."_high"))'
    } else {
      PlotWith="candlesticks whiskerbars"
      PlotUsing='(column(word(FieldNames,fld)."_low")) : (column(word(FieldNames,fld)."_min")) : (column(word(FieldNames,fld)."_max")) : (column(word(FieldNames,fld)."_high"))'
    }
  }
  PlotCmd="plot for [File=1:NumFiles] for [fld=1:NumFields] for [f=fb_min:fb_max] word(PlotFile,File) using ((("
  PlotCmd=PlotCmd.my_x_axis.") == 0 ? 0 : f==0 ? ("
  PlotCmd=PlotCmd.my_x_axis.") : nt - ("
  PlotCmd=PlotCmd.my_x_axis."))+(((File-1)*NumFields+fld-1)*NumFB+f-fb_min)*XF1+XF2) : "
  PlotCmd=PlotCmd.PlotUsing." with "
  PlotCmd=PlotCmd.PlotWith.' title ( SaveFile == 2 ? word(PlotFile,File)." " : "").fb_prefix[f+1].word(FieldNames,fld)'
  #PlotCmd=PlotCmd.PlotWith.' title "ΔT = ".word("12 14 16 20 24 28 32",File)'
  eval PlotCmd
}
EOFMark
}

if [[ "$*" == "" ]];
then
  echo "$0"
  echo "Plot summary data."
  echo "Precede with optional modifiers (key=value):"
  echo "ti     Initial timeslice"
  echo "tf     Final timeslice"
  echo "key    location for key (default: top right)"
  echo "fields Names of fields to display (default: cosh)"
  echo "ref    y-value for reference line"
  echo "err    error for reference line"
  echo "reftext text for reference line"
  echo "log    1 to plot y on log scale, -1 to negate before plotting"
  echo "title  Title for the plot"
  echo "offset X-axis offset between series, 0 to disable (default: 0.05)"
  echo "save   \"1\" to save plots to auto-generated filenames, otherwise name of pdf"
  echo "nt     Plot backward propagating wave as well (or backward only if nt<0)"
  echo "x      definition of field for x-axis. Normally 'column(1)', but try, say 'column(1) / 12'"
  echo "whisker non-zero to plot as box (1 sigma) + error bars (min/max)"
  exit 2
fi

if [[ "$nt" != "" && "$tf" == "" ]]; then tf=$(( (nt < 0 ? -nt : nt)/2 )); fi
if [[ "${title+z}" == "z" ]]; then do_title=1; else do_title=0; fi

# If we specify a filename, then we want a single plot
unset SaveFileName
if [[ "$save" == "" ]]
then
  SaveFile=0
elif [[ "$save" == "1" ]]
then
  SaveFile=1
else
  # A single Filename has been specified - put all the plots in the one file
  SaveFile=2
  SaveFileName="${save// /_}"
  PlotFile="$*"
  PlotFunction
  exit
fi

# Loop through all the files on the command-line performing plots
for PlotFile
do
  PlotPathSplit "$PlotFile"
  if [[ "$mmplotfile_ext" == "txt" ]]
  then #Silently skip non-text files
    #log_limit=1
    #if [[ "$mmplotfile_type" == "corr" ]]; then log_limit=2; fi
    #for(( do_log=0 ; do_log < log_limit ; do_log = do_log + 1 ))
    #do
    PlotFunction
  fi
done
