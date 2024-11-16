#!/usr/bin/env bash

############################################################

# Plot combined (cross ensemble) ZV

############################################################

. PlotCommon.sh NoPreamble

#set -x
set -e

############################################################

# Initialise variables from user

############################################################

if ! [ -d C1 ]; then
  echo "Ensemble C1 doesn't exist. Change directory?"
  exit 2
fi
OutDir=Global/Renorm
ZVDir=${ZVDir:-ZV}
ZVPlotDir=${ZVPlotDir:-${ZVDir}Plot}
ZVEFit=E_For_ZV.txt
ZVFit=ZV.txt

############################################################

# Plot combined (cross ensemble) ZV

############################################################

function PlotZV()
{
gnuplot <<-EOFMark

EnsGroup="$1"
dTMin=$2
Quarks="l $3"
Ensembles="$4"
Seeds="$5"
Spec="$Spec"
Suffix="$Suffix"
Save="${Save}"
dTMax=$dTMax

Field='corr'
FieldL=Field.'_low'
FieldH=Field.'_high'
NumEns=words(Ensembles)
Gap=0.2
OffsetL=(NumEns-1)*0.5*Gap

GetSymb(iEns)=iEns==1 ? 10 : iEns==NumEns ? 8 : 12
InFile(iEns,dT,Q)=word(Ensembles,iEns).'/Renorm/ZV/'.Spec.'p2/ZV_'.Q.'_dt_'.dT.'_'.Suffix.'.fold.'.word(Seeds,iEns).'.txt'
set term pdfcairo dashed size 12in,8in
set key bottom center maxrows NumEns
set pointintervalbox 0 # disables the white space around points in gnuplot 5.4.6 onwards

set linetype 1 ps 0.66 pt 6 lc rgb 0xFA60E4 #pink
set linetype 2 ps 0.66 pt 4 lc rgb 'dark-violet' # 0x9400D3
set linetype 3 ps 0.66 pt 14 lc rgb 0xE36C09 #orange
set linetype 4 ps 0.66 pt 10 lc rgb 0x00B050 #green
set linetype 5 ps 0.66 pt 12 lc rgb 0x0070C0 #blue
set linetype 6 ps 0.66 pt 8 lc rgb 0x003860 #dark navy
set linetype 7 ps 0.66 pt 7 lc rgb 0xFFC000 #yellow
set linetype 8 ps 0.66 pt 7 lc rgb 0xC00000 #red

do for [iQ=1:words(Quarks)] {
Quark=word(Quarks,iQ)
set output Save.'ZV_'.EnsGroup.'_'.Quark.'_spec_'.Spec.'_'.Suffix.'.pdf'

#test
plot for [iEns=0:NumEns-1] for [dT=dTMin:dTMax:4] InFile(iEns+1,dT,Quark) \
  using ( column('t') < 2 || column('t') > dT-2 ? NaN : \
              column('t')-dT*0.5-OffsetL+Gap*iEns+(dT-22)*Gap*0.025) \
        :(column('t') < 2 || column('t') > dT-2 ? NaN : column(Field)) \
        :(column('t') < 2 || column('t') > dT-2 ? NaN : column(FieldL)) \
        :(column('t') < 2 || column('t') > dT-2 ? NaN : column(FieldH)) \
        with yerrorbars \
        lt (dT-8)/4 pt GetSymb(iEns+1) notitle, \
  for [iEns=0:NumEns-1] NaN with points lc 0 pt GetSymb(iEns+1) title word(Ensembles,iEns+1), \
  for [dT=dTMin:dTMax:4] NaN with lines lt (dT-8)/4 title 'ΔT='.dT
set output
}
EOFMark
}

dTMax=24
Save='Global/Plot/'
mkdir -p "$Save"
for Suffix in g5P_g5P g5P_g5W g5W_g5P g5W_g5W; do
for Spec in s l; do
  PlotZV C 12 h6413 'C1 C2'    '2263212701 694228835'
  PlotZV M 16 h447  'M1 M2 M3' '1835672416 2475361470 3301941204'
done
done
