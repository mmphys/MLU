#!/usr/bin/env gnuplot -c

Normalise=1
C2TimeRev=1

if( ARGC < 2 ) {
  print ARG0." Corr1 Corr2 [Output_file [Title [Label1 [Label2 [Negate [field=corr]]]]]"
  print "Compare two correlators, second is time reversed"
  print "Normalise correlators separately on each timeslice so C1(t)=1"
  print "Output [optional] is .pdf to save"
  exit 1
}

NumFiles=2
array Files[NumFiles]
Files[1]=( ARG1 ne "" ) ? ARG1 : 'quark_l_h447_gT_dt_20_p2_0_ps2_1_g5W_g5P.fold.1835672416.txt'
Files[2]=( ARG2 ne "" ) ? ARG2 : 'quark_h447_l_gT_dt_20_p2_1_ps2_0_g5P_g5W.fold.1835672416.txt'

array MyColour[NumFiles]
MyColour[1]=3
MyColour[2]=2

# Define a function to extract deltaT from filename
GetDTSub(s)=strstrt(s,"_dt_") ? substr(s,strstrt(s,"_dt_")+4,strlen(s)) : "0_"
GetDT(s)=substr(GetDTSub(s),1,strstrt(GetDTSub(s),"_")-1)+0

if( C2TimeRev ) {
  DeltaT=GetDT(Files[1])
  DeltaT2=GetDT(Files[2])
  if( DeltaT != DeltaT2 ) {
    print "Delta T ".DeltaT." != Delta T ".DeltaT2
    exit 1
  }
  if( DeltaT == 0 ) {
  DeltaT = 64
  } else {
    if( DeltaT < 32 ) {
      set xrange[-0.5:DeltaT+0.5]
    } else {
      set xrange[0.5:DeltaT-0.5]
    }
  }
}

if( !Normalise ) {
  set logscale y
  set key bottom left
}

# Save to .PDF if third argument specified
if( ARGC >= 3 && ARG3 ne "" ) {
  set term pdf
  set output ARG3.".pdf"
  set pointsize 0.5
}

if( ARGC >=4 && ARG4 ne "" ) {
  set title ARG4
}

array Labels[NumFiles]
Labels[1]=( ARG5 ne "" ) ? ARG5 : Files[1]
Labels[2]=( ARG6 ne "" ) ? ARG6 : Files[2]
IsEnhancedA=( ARG5 ne "" ) ? "enhanced" : "noenhanced"
IsEnhancedB=( ARG6 ne "" ) ? "enhanced" : "noenhanced"

Negate=( ARG7 ne "" ) ? "-" : ""

Field=( ARG8 ne "" ) ? ARG8 : 'corr'
FieldNum=0+system("awk '!/^#/ {for(i=1;i<NF;++i){if($i==\"".Field."\"){print i;exit}}}' ".Files[1])
#FieldsPerFile=0+system("awk '!/^#/ {print NF;exit}' ".Files[1])
#print Field."=".sprintf("%d", FieldNum)." of ".sprintf("%d", FieldsPerFile)

# Set up Awk command to read in the fields we're interested in only
AwkPart1="awk '\\''!/(^#|^t)/ {if( $1 >= 0 && $1 <= ".DeltaT." ) print "
AwkPart2="$1"
do for [i = -1:1] {
  AwkPart2=AwkPart2.",$".sprintf("%d",FieldNum+i) }
AwkPart2=AwkPart2."}'\\'' "

Command="< exec bash -c 'join <(".AwkPart1.AwkPart2.Files[1].") <("
if( C2TimeRev ) {
  Command=Command.AwkPart1.sprintf("%d",DeltaT)."-".AwkPart2.Files[2]."|sort -n"
} else {
  Command=Command.AwkPart1.AwkPart2.Files[2]
}
Command=Command.")'"
#print Command

sNorm=Normalise ? "/column(3)" : ""
UsingC1="(column(3)".sNorm."):(column(2)".sNorm."):(column(4)".sNorm.")"
UsingC2="(".Negate."column(6)".sNorm."):(".Negate."column(5)".sNorm."):(".Negate."column(7)".sNorm.")"

#set yrange [-4:6]

plot Command \
     using 1:@UsingC1 with yerrorbars title Labels[1] @IsEnhancedA lines MyColour[1], \
  '' using ($1+0.1):@UsingC2 with yerrorbars title Labels[2] @IsEnhancedB  lines MyColour[2]
