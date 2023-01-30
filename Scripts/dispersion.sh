#!/bin/sh

if [[ "$*" != "" ]]
then
  echo 'Parameters unexpected'
else
if [[ "$chi" == "" ]]; then do_chi=0; else do_chi=1; fi
if [[ "$timin" == "" ]]; then do_timin=0; else do_timin=1; fi
if [[ "$timax" == "" ]]; then do_timax=0; else do_timax=1; fi
if [[ "$tfmin" == "" ]]; then do_tfmin=0; else do_tfmin=1; fi
if [[ "$tfmax" == "" ]]; then do_tfmax=0; else do_tfmax=1; fi

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark

L=${L:-24}
q1="${q1:-s}"
q2="${q2:-l}"
#seed="${seed:-4147798751}"
seed="${seed:-1835672416}"
#fit="${fit:-corr.g5_gT5}"
fit="${fit:-*}"
do_chi=${do_chi}
chi_max=${chi:-0}
do_timin=${do_timin}
timin=${timin:-0}
do_timax=${do_timax}
timax=${timax:-0}
do_tfmin=${do_tfmin}
tfmin=${tfmin:-0}
do_tfmax=${do_tfmax}
tfmax=${tfmax:-0}
num=${num:-15}
RefVal=${ref:--777}
RefText="E_0=${ref}"
if( "${reftext}" ne "" ) { RefText=RefText." ${reftext}" }
DoOldMomenta=0${oldmom+1}

OutFileName=q1."_".q2
if( fit ne "*" && fit ne "" ) { OutFileName=OutFileName.'.'.fit }
OutFileName=OutFileName.".disp.".seed.".pdf"

Momenta="0 1 2 3 4 5 6"
MomentumSuffix="2"
if( DoOldMomenta ) {
  Momenta="0_0_0 1_0_0 1_1_0 1_1_1 2_0_0 2_1_0 2_1_1"
  MomentumSuffix=""
}
FilenamePrefix=q1."_".q2."_p".MomentumSuffix."_"
FilenameSuffix=".".fit.".params_sort.".seed.".txt"
filename(n)=system("ls ".FilenamePrefix.word(Momenta,n+1).FilenameSuffix)

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

#Work out best fit of E0 in range
TestClause=''
if( do_timin ) { TestClause=TestClause.'column("ti") < TIMin ? 1/0 : ' }
if( do_timax ) { TestClause=TestClause.'column("ti") > TIMax ? 1/0 : ' }
if( do_tfmin ) { TestClause=TestClause.'column("tf") < tfmin ? 1/0 : ' }
if( do_tfmax ) { TestClause=TestClause.'column("tf") > tfmax ? 1/0 : ' }
if( do_chi   ) { TestClause=TestClause.'column("ChiSqPerDof") > chi_max ? 1/0 : ' }
Condition='count > 0 ? 1/0 : '.TestClause
count=0
stats filename(0) using (@Condition (count=count + 1, column("E0"))) nooutput
E0=STATS_min
count=0
stats filename(0) using (@Condition (count=count + 1, column("ti"))) nooutput
WhichTI=int(STATS_min)
count=0
stats filename(0) using (@Condition (count=count + 1, column("tf"))) nooutput
WhichTF=int(STATS_min)
#print 'STATS_min='.gprintf("%g",STATS_min).' STATS_max='.gprintf("%g",STATS_max)

# These are the (components of p)^2 terms in continuum dispersion relation (excluding 0)
array PComponent[2]
do for [i=1:|PComponent|] {
  PComponent[i]=i * 2 * pi / L
  PComponent[i]=PComponent[i] * PComponent[i]
}

# This is p dot p (sum of p^2 terms) in continuum dispersion relation (excluding 0)
array PSquared[6]
PSquared[1]=PComponent[1]
PSquared[2]=PComponent[1] * 2
PSquared[3]=PComponent[1] * 3
PSquared[4]=PComponent[2] + PComponent[1]
PSquared[5]=PComponent[2] + PComponent[1] * 2
PSquared[6]=PComponent[2] + PComponent[1] * 3

# This is the sin^2 term in Lattice dispersion relation (excluding 0)
array SinTerm[2]
do for [i=1:|SinTerm|] {
  SinTerm[i]=sin( i * pi / L )
  SinTerm[i]=SinTerm[i] * SinTerm[i]
}

# This is the sum of the sin^2 terms in Lattice dispersion relation (excluding 0)
array SumSinTerm[6]
SumSinTerm[1]=SinTerm[1]
SumSinTerm[2]=SinTerm[1] * 2
SumSinTerm[3]=SinTerm[1] * 3
SumSinTerm[4]=SinTerm[2]
SumSinTerm[5]=SinTerm[2] + SinTerm[1]
SumSinTerm[6]=SinTerm[2] + SinTerm[1] * 2

MyField="E0"
#Dispersion relation from N^2 (squared momentum in lattice units)
disp(NSquared)=sqrt( E0 * E0 + ( NSquared < 1 ? 0 : PSquared[NSquared >= |PSquared| ? |PSquared| : NSquared] ) )
#Lattice dispersion relation
#Delta=1/L
SinhSqE0=sinh(E0/2)
SinhSqE0=SinhSqE0 * SinhSqE0
LatDisp(NSquared)=NSquared < 1 ? E0 : 2*asinh(sqrt(SinhSqE0 + SumSinTerm[NSquared >= |SumSinTerm| ? |SumSinTerm| : NSquared]))

#sChiDescr='{/Times:Italic χ}^2 per d.o.f.'
sChiDescr='p-value'

set term pdf
set output OutFileName

#set xrange[${xrange:=-0.2:4.2}]
MyTitle="a E_{eff} vs k^2 (".q1."-".q2.")"
if( do_chi ) { MyTitle=MyTitle." (".sChiDescr." ≤ ".sprintf("%g",chi_max).")" }
set title MyTitle
set ylabel 'a E_{eff}'
#set xlabel 'k^2 (sorted by '.sChiDescr.", ".num." best fits)"
set xlabel 'n^2 (sorted by '.sChiDescr.")"
set xtics 1
set key top left reverse Left
set pointsize 0.5

#plot for [i=0:4] filename(i) index "[tf=".tf."]" every ::4::4 using (i):(\$3*\$3):(column(-2)) linecolor variable title "k^2=".i, \
    disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2"

#MyOffset='i + (idx == last_idx ? count : (last_idx=idx, idx==0 ? (count=0,count) : count)) * MyScale'
MyOffset='i + (last_idx == i ? count : (last_idx=i,count=0,count)) * MyScale'
Condition='count >= num ? 1/0 : '.TestClause
#YField='column(MyField)*column(MyField)'
YField='column(MyField)'
#YFieldLow='column(MyField."_low")*column(MyField."_low")'
YFieldLow='column(MyField."_low")'
#YFieldHigh='column(MyField."_high")*column(MyField."_high")'
YFieldHigh='column(MyField."_high")'
MySeries='column("ti") - TIMin + 1'
WithLabels='with labels font "Arial,6" offset char 0,'

#set label "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2" at graph 1, graph 0 font "Arial,12" front textcolor "grey40" offset character -1.5, character 1.5 right
#set label "E_0=0.99656(95), ArXiv:1812.08791" at graph 1, graph 0 font "Arial,8" front textcolor "grey40" offset character -1, character 1 right
if( RefVal == -777 ) { RefText='E_0='.gprintf("%g",E0).' from fit at ti='.WhichTI.', tf='.WhichTF }
set label RefText at graph 1, graph 0 font "Arial,8" front textcolor "grey40" offset character -1, character 1 right
set xrange [0:4.99]
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
