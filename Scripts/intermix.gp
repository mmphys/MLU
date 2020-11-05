#!/usr/bin/env gnuplot

# Interactive mixing
# run from xterm

# Environment variables (add to command-line before calling script)
FieldName="`echo ${field:-exp}`"
Src="`echo ${src:-h3}`"
Snk="`echo ${snk:-h0}`"
phi=`echo ${phi:--90}`
theta=`echo ${theta:--90}`
Normalisation="`echo ${norm:-A1}`"
RefVal=`echo ${ref:-0.19725}`
RefErr=`echo ${referr:-0.00023}`
RefText="`echo ${reftext}`"
PlotType=`echo ${type:-3}`
OtherEnd="`echo ${other:-P}`"

NumDT=4
array DT[NumDT] = [ 12, 14, 16, 20 ]

if( PlotType == 3 ) {
  PlotFile="SinkSourceName"
  TitleFunc="SinkSourceTitle"
} else {
  if( PlotType == 1 ) {
    PlotFile="SourceName"
    TitleFunc="SourceTitle"
  } else {
    if( PlotType == 2 ) {
      PlotFile="SinkName"
      TitleFunc="SinkTitle"
    } else {
      exit error "Environment variable \"type\" invalid"
    }
  }
}

set terminal x11
set xrange [0:1]
set xtics 0.1
set yrange [0.17:0.22]
#set yrange [0.195:0.2]
#set zeroaxis
#plot (x/a)**2, sin(x), 1/x


# Function to give me the filename I need
SinkSourceName(dt,phi,theta)="quark_".Snk."_".Src."_gT_dt_".sprintf("%d",dt)."_p_0_0_0.corr_8_17.corr_8_17.N_".Normalisation.".phi_".sprintf("%d",phi)."_theta_".sprintf("%d",theta)."_g5Pg5W_g5Pg5W.fold.1835672416.txt"
SourceName(dt,phi,theta)="quark_".Snk."_".Src."_gT_dt_".sprintf("%d",dt)."_p_0_0_0.corr_8_17.N_".Normalisation.".theta_".sprintf("%d",theta)."_g5".OtherEnd."_g5Pg5W.fold.1835672416.txt"
SinkName(dt,phi,theta)="quark_".Snk."_".Src."_gT_dt_".sprintf("%d",dt)."_p_0_0_0.corr_8_17.N_".Normalisation.".phi_".sprintf("%d",phi)."_g5Pg5W_g5".OtherEnd.".fold.1835672416.txt"

# Function to give me the title I need
SinkSourceTitle(phi,theta)="Sink and source mixed: phi=".sprintf("%d",phi).", theta=".sprintf("%d",theta)
SinkTitle(phi,theta)="Sink mixed-".OtherEnd." source: phi=".sprintf("%d",phi)
SourceTitle(phi,theta)=OtherEnd." sink-source mixed: theta=".sprintf("%d",theta)

# Setup names of lower and upper limits
FieldLow=FieldName."_low"
FieldHigh=FieldName."_high"


if( RefVal != -777 ) {
  RefErrString=""
  if( RefErr != -777 ) {
    set object 1 rect from graph 0, first RefVal - RefErr to graph 1, first RefVal + RefErr fs solid 0.05 noborder fc rgb "gray10" behind
    RefErrString=sprintf(" (%g)", RefErr)
  }
  set arrow from graph 0, first RefVal to graph 1, first RefVal nohead front lc rgb "gray40" lw 0.25  dashtype "-"
  if( RefText eq "" ) { RefText="Ref: ".sprintf("%g", RefVal).RefErrString }
  set label 2 RefText at screen 0, screen 0 font "Arial,8" front textcolor "grey40" offset character 0.5, 0.5
}


IncFine=1
IncCoarse=10

Finished=0
while( !Finished ) {
  set title Snk." <- ".Src.": ".@TitleFunc(phi,theta)

  plot for [dt=1:NumDT] @PlotFile(DT[dt],phi,theta) \
    using (column("t")/DT[dt]):(column(FieldName)):(column(FieldLow)):(column(FieldHigh)) with yerrorbars title "dt=".sprintf("%d",DT[dt])


  pause mouse keypress
  #if (exists("MOUSE_BUTTON")) { print "Mouse button clicked" } else {
    if (exists("MOUSE_CHAR")) {
      #print "Keypress \"".MOUSE_CHAR."\""
      if( MOUSE_CHAR eq "x" || MOUSE_CHAR eq "X" ) { print " -> Exiting now"; Finished=1 }
      if( MOUSE_CHAR eq "a" || MOUSE_CHAR eq "w" ) { phi=phi+IncFine }
      if( MOUSE_CHAR eq "A" || MOUSE_CHAR eq "W" ) { phi=phi+IncCoarse }
      if( MOUSE_CHAR eq "z" || MOUSE_CHAR eq "s" ) { phi=phi-IncFine }
      if( MOUSE_CHAR eq "Z" || MOUSE_CHAR eq "S" ) { phi=phi-IncCoarse }
      if( MOUSE_CHAR eq "d" || MOUSE_CHAR eq "w" ) { theta=theta+IncFine }
      if( MOUSE_CHAR eq "D" || MOUSE_CHAR eq "W" ) { theta=theta+IncCoarse }
      if( MOUSE_CHAR eq "c" || MOUSE_CHAR eq "s" ) { theta=theta-IncFine }
      if( MOUSE_CHAR eq "C" || MOUSE_CHAR eq "S" ) { theta=theta-IncCoarse }
      if( phi   < -90 ) { phi   = -90 }
      if( phi   >  90 ) { phi   =  90 }
      if( theta < -90 ) { theta = -90 }
      if( theta >  90 ) { theta =  90 }
    }
  #}
}
