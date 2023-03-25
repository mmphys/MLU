#!/usr/bin/env bash

# Optional environment variables
# dir: subdirectory, e.g."frozen"

# Plot ratios
. PlotCommon.sh
PCMakeScriptOutDir

#set -x

# Computed from input

# Optional environment variables
# nodt:  Set to anything to disable plots individual DeltaT all point-wall
# nopw:  Set to anything to disable plots individual point-wall all DeltaT
#DoDeltaT=$((1-0${nodt+1}))
#DoPW=$((1-0${nopw+1}))
Old=$((0${dir:+1}))
series="${series-disp old priorPW betterPW priorP betterP}"
yRange=${yRange:-0.65:1.3}
yRangeRaw=${yRangeRaw:-0.5:1.5}
MaxX=${MaxX:-2.1}

###################################################
# Make a plot of all four form factors - old version - various DeltaT
###################################################

function PlotFunction()
{
  for yAxis in fPlus f0 fPerp fPar
  do
gnuplot <<-EOFMark
Dir="$PlotData/ffs/$1/Summary/"
Meson="$2"
MesonSave="$3"
MyTitle="$4"
RatioNum="$5"
xAxis="${xAxis:-qSq}"
yAxis="${yAxis:-fPlus}"
aInv=${EnsembleaInv:-1.}
xScale=aInv * aInv
MaxPSq=${MaxPSq:-4}

# fPar and fPerp each require different scaling
AdjustYAxis=( yAxis eq "fPerp" || yAxis eq "fPar" )
yAxisLabel=( !AdjustYAxis ) ? "" : ( yAxis eq "fPerp" ) ? " * a^{0.5}" : " * a^{-0.5}"
yScale = ( !AdjustYAxis ) ? 1. : sqrt(aInv)
if( yAxis eq "fPerp" ) { yScale = 1. / yScale }

xAxisH=xAxis."_high"
xAxisL=xAxis."_low"
yAxisH=yAxis."_high"
yAxisL=yAxis."_low"

set title MyTitle.", $Ensemble (R".RatioNum.")"
if( xAxis eq  "qSq" ) { xAxisLabel="q^2" } else { xAxisLabel=xAxis }
if( aInv != 1 ) { xAxisLabel=xAxisLabel." / GeV^2" }
set xlabel xAxisLabel
set ylabel yAxis.yAxisLabel
set key bottom center maxrows 3

f(dt,p)="F".RatioNum."_".Meson."_dt_".dt."_p2_".p.".corr.g5P_g5P.params.1835672416.txt"
#print "f(24,0)=".f(24,0)

DD=0.00025
DD=0.0025

set term pdfcairo font "Arial,12" size 7 in, 3 in
set output "F".RatioNum."_".MesonSave."_".yAxis.".corr.g5P_g5P.pdf"
set pointsize 0.5

if( yAxis eq "fPlus" || yAxis eq "fPerp" ) { pMin=1 } else { pMin=0 }

plot for [p=pMin:MaxPSq] for [dt=24:32:4] Dir.f(dt,p) \
    using (column(xAxis)*xScale+(dt-28)*DD):(column(yAxis)*yScale) \
    :(column(xAxisL)*xScale+(dt-28)*DD):(column(xAxisH)*xScale+(dt-28)*DD) \
    :(column(yAxisL)*yScale):(column(yAxisH)*yScale) with xyerrorbars title "ΔT=".dt.", p^2=".p

EOFMark
  done
}

###################################################
# Make a plot of all four form factors
###################################################

function PlotFuncNew()
{
  for yAxis in fPlus f0 fPerp fPar
  do
gnuplot <<-EOFMark
Dir="$1/"
Meson="$2"
MesonSave="$3"
MyTitle="$4"
RatioNum="$5"
xAxis="${xAxis:-qSq}"
yAxis="${yAxis:-fPlus}"
aInv=${EnsembleaInv:-1.}
xScale=aInv * aInv
MaxPSq=${MaxPSq:-4}
MaxX="$MaxX"
yRange="$yRange"
yRangeRaw="$yRangeRaw"

# fPar and fPerp each require different scaling
AdjustYAxis=( yAxis eq "fPerp" || yAxis eq "fPar" )
yAxisLabel=( !AdjustYAxis ) ? "" : ( yAxis eq "fPerp" ) ? " * a^{0.5}" : " * a^{-0.5}"
yScale = ( !AdjustYAxis ) ? 1. : sqrt(aInv)
if( yAxis eq "fPerp" ) { yScale = 1. / yScale }

xAxisH=xAxis."_high"
xAxisL=xAxis."_low"
yAxisH=yAxis."_high"
yAxisL=yAxis."_low"

set title MyTitle.", $Ensemble (R".RatioNum.")"
if( xAxis eq  "qSq" ) { xAxisLabel="q^2" } else { xAxisLabel=xAxis }
if( aInv != 1 ) { xAxisLabel=xAxisLabel." / GeV^2" }
set xlabel xAxisLabel
set ylabel yAxis.yAxisLabel
set key top left maxrows 3

ModelSuffix=".g5P_g5W.model.1835672416"
f(p)="F".RatioNum."_".Meson."_p2_".p.ModelSuffix.".txt"
#print "f(0)=".f(0)

# Next three lines could be streamlined - see plot command
# Hangover from plotting multiple DeltaT on same plot
DD=0.00025
DD=0.0025
dt=28

set term pdfcairo font "Arial,12" size 7 in, 3 in
set output "F".RatioNum."_".MesonSave."_".yAxis.ModelSuffix.".pdf"
set pointsize 0.5

if( yAxis eq "fPlus" || yAxis eq "fPerp" ) { pMin=1 } else { pMin=0 }

if( MaxX ne "" ) { eval 'set xrange[*:'.MaxX.']' }
if(  AdjustYAxis && yRangeRaw ne "" ) { eval 'set yrange['.yRangeRaw.']' }
if( !AdjustYAxis && yRange    ne "" ) { eval 'set yrange['.yRange.']' }

plot for [p=MaxPSq:pMin:-1] Dir.f(p) \
    using (column(xAxis)*xScale+(dt-28)*DD):(column(yAxis)*yScale) \
    :(column(xAxisL)*xScale+(dt-28)*DD):(column(xAxisH)*xScale+(dt-28)*DD) \
    :(column(yAxisL)*yScale):(column(yAxisH)*yScale) with xyerrorbars \
    linestyle p+1 title "n^2=".p

EOFMark
  done
}

###################################################
# Main
###################################################

#Spec=(sp2 lp2 lp2)
Spec=(sp2)
Meson=(l_h$Heavy s_h$Heavy l_h$Heavy)
MesonSave=(K_Ds$Heavy K_D$Heavy pi_D$Heavy)
Title=("D_s ⟹ K" "D ⟹ K" "D ⟹ π")

if (( Old ))
then
  mkdir -p $dir
  cd $dir
for (( i=0; i < ${#Spec[@]}; ++i ))
do
  echo $i/${#Spec[@]} PlotFunction "${Spec[i]}" "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 1
                      PlotFunction "${Spec[i]}" "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 1
  if (( i==0 )); then PlotFunction "${Spec[i]}" "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 3; fi
done
fi

for OutSub in $series; do
  mkdir -p $OutSub
  cd $OutSub
  i=0
  InPrefix=../.. # $HOME/NoSync/$Ensemble
  MELFit=$InPrefix/MELFit
  SpecDir="3sm_${Spec[i]}"
  FitFile=$MELFit/Fit_${Spec[i]}_$OutSub.txt
  Cmd="CRatio --type f,$L --efit $FitFile --i3 $MELFit/$SpecDir/ -o $SpecDir/ '*.$OutSub.*.h5'"
  LogFile="FFS_$SpecDir.log"
  echo "$Cmd"
  echo "$Cmd"   > $LogFile
  eval "$Cmd" &>> $LogFile
  echo PlotFuncNew $SpecDir "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 3
       PlotFuncNew $SpecDir "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 3
  cd ..
done

#series=all yRange=0.65:1.3 yRangeRaw=0.5:1.5 MaxX=2.1 PlotFormFactor.sh
