#!/usr/bin/env bash

# Optional environment variables
# dir: subdirectory, e.g."frozen"

# Plot ratios
. PlotCommon.sh
mkdir -p All/FormFactor
cd All/FormFactor

#set -x

# Computed from input
OutSub=${dir:-final}

# Optional environment variables
# nodt:  Set to anything to disable plots individual DeltaT all point-wall
# nopw:  Set to anything to disable plots individual point-wall all DeltaT
#DoDeltaT=$((1-0${nodt+1}))
#DoPW=$((1-0${nopw+1}))

mkdir -p $OutSub; cd $OutSub

###################################################
# Make a plot of all four form factors
###################################################

function PlotFunction()
{
  for yAxis in fPlus f0 fPerp fPar
  do
gnuplot <<-EOFMark
Dir="/analyse/ffs/$1/Summary/"
Meson="$2"
MesonSave="$3"
MyTitle="$4"
RatioNum="$5"
DoDisp=0$6
DeltaT=$DeltaT
xAxis="${xAxis:-qSq}"
yAxis="${yAxis:-fPlus}"

if( DoDisp ) { xAxis=xAxis."Lat" }

PlotData="/Volumes/QCD/tursa/semilep/data"
Ensembles="F1M M1 M2 M3 C1 C2"
NumEnsembles=words(Ensembles)
HeavyMasses="385 447 447 447 6413 6413"

array EnsembleaInv[6]
EnsembleaInv[1]=2.708
EnsembleaInv[2]=2.3833
EnsembleaInv[3]=2.3833
EnsembleaInv[4]=2.3833
EnsembleaInv[5]=1.7848
EnsembleaInv[6]=1.7848

array EnsembleaInvSq[6]
do for [idx=1:|EnsembleaInv|] { EnsembleaInvSq[idx] = EnsembleaInv[idx] * EnsembleaInv[idx] }

array pMax[6]
pMax[1]=6
pMax[2]=4
pMax[3]=4
pMax[4]=4
pMax[5]=4
pMax[6]=4

# fPar and fPerp each require different scaling
array EnsYScale[6]
AdjustYAxis=( yAxis eq "fPerp" || yAxis eq "fPar" )
yAxisLabel=( !AdjustYAxis ) ? "" : ( yAxis eq "fPerp" ) ? " * a^{0.5}" : " * a^{-0.5}"
do for [idx=1:|EnsYScale|] {
    EnsYScale[idx] = ( !AdjustYAxis ) ? 1. : sqrt(EnsembleaInv[idx])
    if( yAxis eq "fPerp" ) { EnsYScale[idx] = 1. / EnsYScale[idx] }
}

xAxisH=xAxis."_high"
xAxisL=xAxis."_low"
yAxisH=yAxis."_high"
yAxisL=yAxis."_low"

# Split Meson into
MesonPos=strstrt(Meson, "*")
if( MesonPos ) {
  MesonL=substr(Meson,1,MesonPos-1)
  MesonR=substr(Meson,MesonPos+1,strlen(Meson))
} else {
  MesonL=Meson
  MesonR=""
}

MyTitle=MyTitle."      (R".RatioNum.")"
MyTitle=MyTitle."      ΔT=".DeltaT
set title MyTitle
if( xAxis eq  "qSq" || xAxis eq "qSqLat" ) { xAxisLabel="q^2" } else { xAxisLabel=xAxis }
xAxisLabel=xAxisLabel." / GeV^2"
set xlabel xAxisLabel
set ylabel yAxis.yAxisLabel
set key bottom right maxrows 3
if( MesonSave eq "K_Ds" ) {
  set xrange [-0.2:2.3]
  if( RatioNum == 3 ) {
    set yrange [0.5:1.2]
  } else {
    set yrange [0.2:1.4]
  }
}

f(dt,p)="F".RatioNum."_".MesonL.(MesonPos ? word(HeavyMasses, Ens) : "").MesonR. \
        "_dt_".dt."_p2_".p.".corr.g5P_g5P.params.1835672416.txt"
#print "f(24,0)=".f(24,0)

DeltaTStyle(dt)=dt/2-8

array EnsembleColour[6]
EnsembleColour[1]=7
EnsembleColour[2]=3
EnsembleColour[3]=1
EnsembleColour[4]=2
EnsembleColour[5]=5
EnsembleColour[6]=4

OutFileName="F".RatioNum."_".MesonSave."_dt_".DeltaT."."
if( DoDisp ) { OutFileName=OutFileName."LatDisp" } else { OutFileName=OutFileName."EFit" }
OutFileName=OutFileName.".".yAxis.".corr.g5P_g5P.pdf"

set term pdfcairo font "Arial,12" size 7 in, 3 in
set output OutFileName
set pointsize 0.5

if( yAxis eq "fPlus" || yAxis eq "fPerp" ) { pMin=1 } else { pMin=0 }

if( DoDisp ) {
  MyUsing=""
  MyWith="yerrorbars"
} else {
  MyUsing=":(column(xAxisL)*EnsembleaInvSq[Ens]):(column(xAxisH)*EnsembleaInvSq[Ens])"
  MyWith="xyerrorbars"
}

plot \
  for [Ens=1:NumEnsembles] for [p=pMin:pMax[Ens]] \
    PlotData."/".word(Ensembles, Ens).Dir.f(DeltaT,p) \
    using (column(xAxis)*EnsembleaInvSq[Ens]) \
      :(column(yAxis) * EnsYScale[Ens]) \
      @MyUsing \
      :(column(yAxisL) * EnsYScale[Ens]) \
      :(column(yAxisH) * EnsYScale[Ens]) \
    with @MyWith \
    linecolor EnsembleColour[Ens] pointtype DeltaTStyle(DeltaT) notitle, \
  for [Ens=1:NumEnsembles] '' every ::::0 \
    using (NaN):(NaN) \
    with lines \
    linecolor EnsembleColour[Ens] title word(Ensembles, Ens)#, \
  for [dt=24:32:4] '' every ::::0 \
    using (NaN):(NaN) \
    with points \
    pointtype DeltaTStyle(dt) title "ΔT=".dt

EOFMark
  done
}

###################################################
# Main
###################################################

Spec=(sp2 lp2 lp2)
#Spec=(sp2)
Meson=('l_h*' 's_h*' 'l_h*')
MesonSave=(K_Ds K_D pi_D)
Title=("D_s ⟹ K" "D ⟹ K" "D ⟹ π")

for Disp in 0 1
do
for DeltaT in {24..32..4}
do
  for (( i=0; i < ${#Spec[@]}; ++i ))
  do
    echo $i/${#Spec[@]} PlotFunction "${Spec[i]}" "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 1 $Disp
                        PlotFunction "${Spec[i]}" "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 1 $Disp
    if (( i==0 )); then PlotFunction "${Spec[i]}" "${Meson[i]}" "${MesonSave[i]}" "${Title[i]}" 3 $Disp; fi
  done
done
done
