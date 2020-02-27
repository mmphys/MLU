#!/bin/sh

if [[ "$*" != "" ]]
then
  echo 'Parameters unexpected'
else
if [[ "$chi" == "" ]]; then do_chi=0; else do_chi=1; fi
if [[ "$timin" == "" ]]; then do_timin=0; else do_timin=1; fi
if [[ "$timax" == "" ]]; then do_timax=0; else do_timax=1; fi

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark

q1="${q1:-h1}"
q2="${q2:-l}"
seed="${seed:-4147798751}"
fit="${fit:-corr.g5_gT5}"
do_chi=${do_chi}
chi_max=${chi:-0}
do_timin=${do_timin}
timin=${timin:-0}
do_timax=${do_timax}
timax=${timax:-0}
num=${num:-15}

MyField="E0"
E0=0.99656
L=24
kFactorSq=2 * pi / L
kFactorSq=kFactorSq * kFactorSq
#Momentum from k_x^2 + k_y^2 + k_z^2
#Momentum(KDotK)=kFactor * sqrt(KDotK)
#Dispersion relation from p
disp2(KDotK)=KDotK >= 5 ? NaN : E0 * E0 + kFactorSq * ( KDotK > 5 ? 5 : int(KDotK) )
disp(KDotK)=sqrt(disp2(KDotK))
#Lattice dispersion relation
#Delta=1/L
SinhSqE0=sinh(E0/2)
SinhSqE0=SinhSqE0 * SinhSqE0
K1Factor=sin(pi/L)
K1Factor=K1Factor*K1Factor
K2Factor=sin(2*pi/L)
K2Factor=K2Factor*K2Factor
SinSum(KDotK)=KDotK < 1 ? 0 : KDotK < 4 ? int( KDotK ) * K1Factor : K2Factor + int(KDotK-4) * K1Factor;
LatDisp(KDotK)=KDotK >= 5 ? NaN : 2*asinh(sqrt(SinhSqE0 + SinSum(KDotK)))
LatDisp2(KDotK)=LatDisp(KDotK)*LatDisp(KDotK)

filename(n)=q1."_".q2."_p_".word("0_0_0 1_0_0 1_1_0 1_1_1 2_0_0",n+1).".".fit.".params_sort.".seed.".txt"

SeriesName(n)="ti=".n

TIMin = 0
TIMax = 0
#set xrange [*:${chi:-*}]
do for [i=0:4] {
  stats filename(i) using "ti" nooutput
  #print "i=".i." STATS_min=".sprintf("%g",STATS_min)." STATS_max=".sprintf("%g",STATS_max)
  # Keep track of overall lowest and highest
  if( i == 0 ) {
    TIMin = STATS_min
    TIMax = STATS_max
  } else {
    if( TIMin > STATS_min ) { TIMin = STATS_min }
    if( TIMax < STATS_max ) { TIMax = STATS_max }
  }
}
if( do_timin ) { TIMin = TIMin < timin ? timin : TIMin }
if( do_timax ) { TIMax = TIMax > timax ? timax : TIMax }
#print "TIMin=".sprintf("%g",TIMin)." TIMax=".sprintf("%g",TIMax)

sChiDescr='{/Times:Italic χ}^2 per d.o.f.'

set term pdf
set output q1."_".q2.".".fit.".disp.4147798751.pdf"

#set xrange[${xrange:=-0.2:4.2}]
MyTitle="a E_{eff} vs k^2 (".q1."-".q2.")"
if( do_chi ) { MyTitle=MyTitle." (".sChiDescr." ≤ ".sprintf("%g",chi_max).")" }
set title MyTitle
set ylabel 'a E_{eff}'
#set xlabel 'k^2 (sorted by '.sChiDescr.", ".num." best fits)"
set xlabel 'k^2 (sorted by '.sChiDescr.")"
set xtics 1
set key top left reverse Left
set pointsize 0.5

#plot for [i=0:4] filename(i) index "[tf=".tf."]" every ::4::4 using (i):(\$3*\$3):(column(-2)) linecolor variable title "k^2=".i, \
    disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2"

#MyOffset='i + (idx == last_idx ? count : (last_idx=idx, idx==0 ? (count=0,count) : count)) * MyScale'
MyOffset='i + (last_idx == i ? count : (last_idx=i,count=0,count)) * MyScale'
Condition='count >= num ? 1/0 : '
if( do_timin ) { Condition=Condition.'column("ti") < TIMin ? 1/0 : ' }
if( do_timax ) { Condition=Condition.'column("ti") > TIMax ? 1/0 : ' }
if( do_chi ) { Condition=Condition.'column("ChiSqPerDof") > chi_max ? 1/0 : ' }
#YField='column(MyField)*column(MyField)'
YField='column(MyField)'
#YFieldLow='column(MyField."_low")*column(MyField."_low")'
YFieldLow='column(MyField."_low")'
#YFieldHigh='column(MyField."_high")*column(MyField."_high")'
YFieldHigh='column(MyField."_high")'
MySeries='column("ti") - TIMin + 1'
WithLabels='with labels font "Arial,6" offset char 0,'

#set label "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2" at graph 1, graph 0 font "Arial,12" front textcolor "grey40" offset character -1.5, character 1.5 right
set label "E_0=0.99656(95), ArXiv:1812.08791" at graph 1, graph 0 font "Arial,8" front textcolor "grey40" offset character -1, character 1 right
set xrange [0:6]
set samples 1000
MyScale=1./num < 0.08 ? 0.08 : 1./num
count = 0
last_idx=-1
plot \
  disp(x) lc rgb "gray20" lw 0.25 dashtype "-" title "Continuum disp", \
  LatDisp(x) lc rgb "#FF1493" lw 0.5 title "Lattice disp", \
  for [i=0:4] filename(i) using \
    (@MyOffset):(@Condition @YField):(@Condition @YFieldLow):(@Condition (count=count+1, @YFieldHigh)):(@MySeries) \
    with yerrorbars notitle lc variable pt 6 pointsize 0.5, \
  for [idx=TIMin:TIMax] filename(0) index 0 using 1:(1/0):(1/0):(1/0) \
    with yerrorbars title SeriesName(idx) lc idx - TIMin + 1 pt 6, \
  for [i=0:4] filename(i) using \
    (@MyOffset):(@Condition (count=count+1, @YFieldHigh)):"tf" @WithLabels 0.5 notitle
#plot \
for [i=0:0] for[idx=0:NBlock[i+1]-1] filename(i) index idx using \
(@MyOffset):(@Condition (count=count+1, @YField)):(@Condition @YFieldLow):(@Condition @YFieldHigh) \
with yerrorbars title columnheader(1) ls idx + 1, \
for [i=1:4] for[idx=0:NBlock[i+1]-1] filename(i) index idx using \
(@MyOffset):(@Condition (count=count+1, @YField)):(@Condition @YFieldLow):(@Condition @YFieldHigh) \
with yerrorbars notitle ls idx + 1, \
for [i=0:4] for[idx=0:NBlock[i+1]-1] filename(i) index idx using \
(@MyOffset):(@Condition (count=count+1, @YFieldHigh)):"tf" @WithLabels 0.5 notitle, \
disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion"

EOFMark

fi
