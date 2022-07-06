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
FileDT="${FileDT}"
my_xrange="${ti:-*}:${tf:-*}"
my_yrange="${yrange:-*:*}"
my_cbrange="${cbrange:--1:1}"
my_ylabel="${ylabel}"
my_xlabel="${xlabel}"
my_key="${key:-top right}"
FieldNames="${fields:-cosh}"
RefVal=${ref:--777}
RefErr="$err" #1 value for relative error, 2 values for low and high
RefText="${reftext-MagicBanana}"
do_log=${log:-0}
SaveFile=${SaveFile:-0}
SaveFileName="$SaveFileName"
SaveLabel="${savelabel-$SaveFileName}"
do_offset=${offset:-0.05}
do_whisker=0${whisker+1} #non-zero if whisker set to anything (including nothing)
do_rel=0${rel}
got_legend=${got_legend}
$PlotFileLegend
PDFSize="${size}"
my_xtics="${xtics}"
my_ytics="${ytics}"
LegEnh="${legenh+enhanced}"
Negate="$negate"
XRevT="$xrevt"

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

#while(words(my_x_axis) < NumFiles) { my_x_axis=my_x_axis." ".word(my_x_axis,words(my_x_axis)) }
array XRevA[NumFiles]
do for [i=1:NumFiles] {
  if( i <= words(XRevT) ) {
    eval 'XRevA[i]='.word(XRevT,i)
  } else {
    XRevA[i]=0
  }
  #print sprintf("XRevA[%d]=",i).XRevA[i]
}
GetXAxis(X,File)=XRevA[File] == 0 ? X : XRevA[File] - X
my_x_axis='GetXAxis('.my_x_axis.',File)'

# Work out whether we have absolute minimum ... and can therefore do box and whiskers
if( do_whisker ) {
do_whisker=(system("awk '! /#/ {print (\$0 ~ /".word(FieldNames,1)."_min/) ? 1 : 0; exit}' ".word(PlotFile,1)) eq "0") ? 0 : 1 }

# Decide on a title and output filename
if( SaveFile == 2 ) {
  MyTitle=my_title
  SaveFileName=SaveFileName.".pdf"
} else {
  # We are only processing one file
  MyTitle="${mmplotfile_name_no_ext}".FBText
  SaveFileName="${mmplotfile_base}.".FileType
  if( FileType ne "cormat" ) { SaveFileName=SaveFileName.FieldNamesShort.LogShort.FBShort }
  SaveFileName=SaveFileName.".${mmplotfile_seed}.pdf"
}

if( SaveFile ) {
  if( PDFSize ne "" ) { PDFSize="size ".PDFSize } else {
    if( NumFiles==1 && FileType eq "cormat" ) { PDFSize="size 7.25,7 fontscale 0.4" }
  }
  eval "set term pdfcairo ".PDFSize
  set pointsize 0.4
  set output SaveFileName
  if( SaveLabel ne "" ) {
    set label 1 SaveLabel noenhanced at screen 1, 0.5 center rotate by -90 \
        font ",8" front textcolor "grey40" offset character -1.3, 0
  }
}

set key @my_key noenhanced
set xrange[@my_xrange]
set yrange[@my_yrange]
if( do_title ) {
  # This title might be set to empty string, meaning - no title
  if( my_title ne "" ) { set title my_title font ', 18' }
} else {
  set title MyTitle noenhanced
}
if( my_xlabel ne "" ) { set xlabel my_xlabel font ', 16' }
if( my_ylabel ne "" ) { set ylabel my_ylabel font ', 16' }

if( RefVal != -777 ) {
  RefErrString=""
  RefXUnits="graph"
  RefXLeft=0
  RefXRight=1
  if( RefErr ne "" ) {
    if( words(RefErr) == 1 ) {
      # 1 word - this is (symmetric) relative error
      RefErrLow=RefVal-RefErr
      RefErrHigh=RefVal+RefErr
      RefErrString=" (".RefErr.")"
    } else {
      if( words(RefErr) == 2 ) {
        # 2 words - these are the absolute values of low and high error bars
        RefErrLow=0.+word(RefErr,1)
        RefErrHigh=0.+word(RefErr,2)
      } else {
        # 4 words - x y (start) x y (end)
        RefErrLow=0.+word(RefErr,2)
        RefErrHigh=0.+word(RefErr,4)
        RefXUnits="first"
        RefXLeft=0.+word(RefErr,1)
        RefXRight=0.+word(RefErr,3)
      }
      RefErrString=sprintf(" (-%g, +%g)", RefVal - RefErrLow, RefErrHigh - RefVal )
    }
    eval "set object 1 rect from ".RefXUnits." RefXLeft, first RefErrLow to ".RefXUnits." RefXRight, first RefErrHigh" \
      ." fs solid 1 noborder fc rgb 0xD0D0D0 behind"
  }
  eval "set arrow from ".RefXUnits." RefXLeft, first RefVal to ".RefXUnits." RefXRight, first RefVal nohead front lc rgb \"gray40\" lw 0.25 dashtype \"-\""
  if( RefText eq "MagicBanana" ) { RefText="Ref: ".sprintf("%g", RefVal).RefErrString }
  #set label 2 RefText at screen 0, screen 0 font ",8" front textcolor "grey40" offset character 0.5, 0.25
}

#AbsMin(y,low,high)=sgn(y) < 0 ? -(high) : low
#AbsMax(y,low,high)=sgn(y) < 0 ? -(low) : high

if( SaveFile != 2 ) {
  # Don't do this for multi-plots, as probably best to stick to colour scheme
  set linetype 1 lc rgb 'blue'
  set linetype 2 lc rgb 'red'
} else {
  # This is a multi-series plot
  set linetype 1 pt 7 lc rgb 0x0070C0 #blue
  set linetype 2 pt 7 lc rgb 0xE36C09 #orange
  set linetype 3 pt 7 lc rgb 0x00B050 #green
  set linetype 4 pt 7 lc rgb 'dark-violet'
  set linetype 5 pt 7
  set linetype 6 pt 7
  set linetype 7 pt 7
  set linetype 8 pt 7
}

# The logging variable becomes a log
if( do_log ) { set logscale y }

# Function to get a field - potentially negated
GetField(Suffix)='('.'(word(Negate,File) eq "-" ? -1 : 1)*column(word(FieldNames,fld)."'.Suffix.'"))'
# Function to get the x-axis for a file
#GetXAxis()=''

# Work out how much to offset each series by
NumFields=words(FieldNames)
NumFB=fb_max - fb_min + 1
XF1=do_offset
XF2=(1 - NumFields*NumFB*NumFiles)/2*XF1

PlotWith="with yerrorbars"
PlotUsing=""
if( NumFiles==1 && FileType eq "cormat" ) {
  eval "set xtics rotate noenhanced ".my_xtics
  eval "set ytics noenhanced ".my_ytics
  set cbrange [@my_cbrange]
  plot PlotFile matrix columnheaders rowheaders with image pixels
} else {
  if( my_xtics ne "" ) { eval "set xtics ".my_xtics }
  if( my_ytics ne "" ) { eval "set ytics ".my_ytics }
#  if( do_log ) {
#    PlotUsing='(abs(column(word(FieldNames,fld)))) : (sgn(column(word(FieldNames,fld)))*column(word(FieldNames,fld).( sgn(column(word(FieldNames,fld))) ? "_low" : "_high" ))) : (sgn(column(word(FieldNames,fld)))*column(word(FieldNames,fld).( sgn(column(word(FieldNames,fld))) ? "_high" : "_low" )))'
#  } else {
    if( do_whisker ) {
      PlotWith="with candlesticks whiskerbars"
      PlotUsing=GetField("_low").' : '.GetField("_min").' : '.GetField("_max").' : '.GetField("_high")
    } else {
      if( do_rel ) {
        if( do_rel == 1 ) {
          PlotUsing='(1):('.GetField("_low").'/'.GetField("").'):('.GetField("_high").'/'.GetField("").')'
        } else {
          PlotWith=""
          PlotUsing='(('.GetField("_high").' - '.GetField("_low").')/(2*'.GetField("").'))'
        }
      } else {
        PlotUsing=GetField("").' : '.GetField("_low").' : '.GetField("_high")
        }
    }
#  }
  #print 'PlotUsing="'.PlotUsing.'"'
  PlotCmd="plot for [File=1:NumFiles] for [fld=1:NumFields] for [f=fb_min:fb_max] word(PlotFile,File) using ((("
  PlotCmd=PlotCmd.my_x_axis.") == 0 ? (f==0 ? 0 : 1/0) : f==0 ? ("
  PlotCmd=PlotCmd.my_x_axis.") : nt - ("
  PlotCmd=PlotCmd.my_x_axis."))+(((File-1)*NumFields+fld-1)*NumFB+f-fb_min)*XF1+XF2) : "
  PlotCmd=PlotCmd.PlotUsing.' '.PlotWith
  PlotCmd=PlotCmd.' title ( (SaveFile == 2 || got_legend) ? PlotFileLegend[File]." " : "").fb_prefix[f+1].(got_legend ? "" : word(FieldNames,fld)) '.LegEnh.' '
  #PlotCmd=PlotCmd.PlotWith.' title "ΔT = ".word("12 14 16 20 24 28 32",File)'
  if( RefVal != -777 && RefText ne "MagicBanana" ) {
    PlotCmd=PlotCmd.", 1/0 lt 0 dashtype 2 title '".RefText."'"
  }
  #print PlotCmd
  eval PlotCmd
  if( SaveFile == 0 ) { pause mouse close }
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
  echo "err    error for reference line e.g. \"0.5\" for +/- 0.5 or low and high values"
  echo "reftext text for reference line"
  echo "log    1 to plot y on log scale, (-1 to negate before plotting - not true any more?)"
  echo "rel    plot relative error 1=error bars, 2=error (does nothing on log scale)"
  echo "title  Title for the plot"
  echo "offset X-axis offset between series, 0 to disable (default: 0.05)"
  echo "save   \"1\" to save plots to auto-generated filenames, otherwise name of pdf"
  echo "savelabel Filename to place on RHS (default=save file name)"
  echo "nt     Plot backward propagating wave as well (or backward only if nt<0)"
  echo "x      definition of field for x-axis. Normally 'column(1)', but try, say 'column(1) / 12'"
  echo "xargs  Arguments for field for x-axis, eg deltaT"
  echo "whisker plot as box (1 sigma) + error bars (min/max)"
  echo "legend string containing legend for each file"
  echo "       e.g. \"'big banana' hulahoop\""
  echo "legenh 1 for enhanced legend"
  echo "xlabel X-axis label"
  echo "ylabel X-axis label"
  echo "size   PDF size command (e.g. \"5.75,1.75 font 'Times-New-Roman,12'\""
  echo "negate One string per file. If any set to '-' then the field value is negated"
  echo "xrevt  One value per file. Any non-zero value means t -> xrevt - t"
  exit 2
fi

if [[ "$nt" != "" && "$tf" == "" ]]; then tf=$(( (nt < 0 ? -nt : nt)/2 )); fi
if [[ "${title+z}" == "z" ]]; then do_title=1; else do_title=0; fi

if [[ "$legend" == "" ]]; then
  got_legend=0
  legend=()
else
  got_legend=1
  eval legend=(${legend})
fi

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
  if (( do_title == 0 )); then title="${SaveFileName##*/}"; fi
  PlotFile="$*"
  for f in $PlotFile; do
    PlotPathSplit "$f"
    FileDT="$FileDT $mmplotfile_dt"
    if !(( got_legend )); then
      legend+=( "${f%.$mmplotfile_seed.$mmplotfile_ext}" )
    fi
  done
  PassGnuPlotArray legend PlotFileLegend
  PlotFunction
  exit
fi

# Loop through all the files on the command-line performing plots
LegendIndex=0
for PlotFile
do
  PlotPathSplit "$PlotFile"
  if [[ "$mmplotfile_ext" == "txt" ]]
  then #Silently skip non-text files
    PlotFileTitle="${PlotFile%.$mmplotfile_seed.$mmplotfile_ext}"
    FileDT="$mmplotfile_dt"
    if (( got_legend )); then
      ThisLegend=( "${legend[$((LegendIndex++))]}" )
    else
      ThisLegend=( "${f%.$mmplotfile_seed.$mmplotfile_ext}" )
    fi
    #log_limit=1
    #if [[ "$mmplotfile_type" == "corr" ]]; then log_limit=2; fi
    #for(( do_log=0 ; do_log < log_limit ; do_log = do_log + 1 ))
    #do
    PassGnuPlotArray ThisLegend PlotFileLegend
    PlotFunction
  fi
done
