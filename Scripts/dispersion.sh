#!/usr/bin/env bash

# Perform a fit for matrix elements given a choice of 2pt fit ranges

############################################################

# Inputs

############################################################

#Ensemble=${Ensemble:-F1M}
. PlotCommon.sh

#set -x
set -e

if [[ "$*" != "" ]]
then
  echo 'Parameters unexpected. Optional environment variables:'
  echo 'q1      Quark 1, default: s'
  echo 'q2      Quark 2, default: l'
  echo 'NumExp  Number of energy levels, default: 2'
  echo 'tdisp   Fit times specifying ground state, e.g. 6_22_7_23 for 2-model fit'
else
  q1="${q1:-s}"
  q2="${q2:-l}"
  OptionNoMass= GetMeson Meson $q1 $q2
  # Parse and validate fit times specifying which fit is ground state in dispersion relation
  DispTI=${tdisp%%_*}
  DispTFLabel=${tdisp#*_}
  if ! [[ -v tdisp && "$tdisp" != "$DispTI" && "$DispTI" && "$DispTFLabel" ]]; then
    DispTI=0
    unset DispTFLabel
  fi

###################################################
# Make a plot of the simulated correlator masses
###################################################

gnuplot <<-EOFMark

L=${L:-24}
MaxPSq=${MaxPSq:-4}
q1="${q1}"
q2="${q2}"
NumExp=${NumExp:-2}
#MLUSeed="${MLUSeed:-4147798751}"
MLUSeed="${MLUSeed:-1835672416}"
#fit="${fit:-corr.g5_gT5}"
fit="${fit:-*}"
do_chi=0${chi+1}
chi_max=${chi:-0}
do_timin=0${timin+1}
timin=${timin:-0}
do_timax=0${timax+1}
timax=${timax:-0}
do_tfmin=0${tfmin+1}
tfmin=${tfmin:-0}
do_tfmax=0${tfmax+1}
tfmax=${tfmax:-0}
num=${num:-15}
RefVal=${ref:--777}
RefText="E_0=${ref}"
if( "${reftext}" ne "" ) { RefText=RefText." ${reftext}" }
DoOldMomenta=0${oldmom+1}
Meson="$Meson"
DispTI=$DispTI
DispTFLabel="$DispTFLabel"

OutFileName=q1."_".q2
if( fit ne "*" && fit ne "" ) { OutFileName=OutFileName.'.'.fit }

Momenta="0 1 2 3 4 5 6"
MomentumSuffix="2"
if( DoOldMomenta ) {
  Momenta="0_0_0 1_0_0 1_1_0 1_1_1 2_0_0 2_1_0 2_1_1"
  MomentumSuffix=""
}
FilenamePrefix=q1."_".q2."_p".MomentumSuffix."_"
FilenameSuffix=".".fit.".params_sort.".MLUSeed.".txt"
filename(n)=system("ls ".FilenamePrefix.word(Momenta,n+1).FilenameSuffix)

SeriesName(n)="ti=".n

TIMin = 0
TIMax = 0
#set xrange [*:${chi:-*}]
do for [i=0:MaxPSq] {
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
# Honour limits on ti min and max
if( do_timin ) { TIMin = TIMin < timin ? timin : TIMin }
if( do_timax ) { TIMax = TIMax > timax ? timax : TIMax }
#print "TIMin=".sprintf("%g",TIMin)." TIMax=".sprintf("%g",TIMax)

#Work out best fit of E0 in range
# Relies on file already being sorted by pvalue (Hotelling) so that first meeting criteria is best
TestClause=''
if( do_timin ) { TestClause=TestClause.'column("ti") < TIMin ? 1/0 : ' }
if( do_timax ) { TestClause=TestClause.'column("ti") > TIMax ? 1/0 : ' }
if( do_tfmin ) { TestClause=TestClause.'column("tf") < tfmin ? 1/0 : ' }
if( do_tfmax ) { TestClause=TestClause.'column("tf") > tfmax ? 1/0 : ' }
if( do_chi   ) { TestClause=TestClause.'column("ChiSqPerDof") > chi_max ? 1/0 : ' }

# Find the first in-range fit
if( DispTFLabel eq "" ) {
  print "*** Begin ignore warnings - finding fit times"
  count=0
  set table
  plot filename(0) using (count > 0 ? 1/0 : @TestClause (count=count + 1, DispTI=int(column("ti"))))\
                        :(count!= 1 ? 1/0 : ( DispTFLabel=stringcolumn("tfLabel"), count=count+1))
  unset table
  print "*** End ignore warnings - finding fit times"
}
#print "Dispersion TI=".DispTI.", TFLabel=".DispTFLabel

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

#Dispersion relation from N^2 (squared momentum in lattice units)
disp(NSquared)=sqrt( E0 * E0 + ( NSquared < 1 ? 0 : PSquared[NSquared >= |PSquared| ? |PSquared| : NSquared] ) )
#Lattice dispersion relation
E0=0
SinhSqE0=0
LatDisp(NSquared)=NSquared < 1 ? E0 : 2*asinh(sqrt(SinhSqE0 + SumSinTerm[NSquared >= |SumSinTerm| ? |SumSinTerm| : NSquared]))

#sChiDescr='{/Times:Italic χ}^2 per d.o.f.'
sChiDescr='p-value'
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
WithLabels='with labels font "Arial,6" rotate noenhanced left offset char 0,'

do for [ExState=0:NumExp - 1] {
MyField="E".ExState
MyFieldLabel="E_{".ExState."}"

#set table
#plot filename(0) using 0:((column("ti")==DispTI && stringcolumn("tfLabel") eq DispTFLabel) ? E0=column(MyField):column(MyField))
#unset table
eval 'stats filename(0) using ((column("ti")==DispTI && stringcolumn("tfLabel") eq DispTFLabel) ? E0=column("'.MyField.'"):1/0) nooutput'
E0=STATS_min
#print 'E0='.gprintf("%g",E0)

#count=0
#stats filename(0) using (@Condition (count=count + 1, column(MyField))) nooutput
#E0=STATS_min
SinhSqE0=sinh(E0/2)
SinhSqE0=SinhSqE0 * SinhSqE0

set term pdf
set output OutFileName.".disp_".MyField.".".MLUSeed.".pdf"

#set xrange[${xrange:=-0.2:4.2}]
#MyTitle="a E_{eff} vs n^2 (".q1."-".q2.")"
MyTitle=Meson." dispersion relation ".MyFieldLabel
if( do_chi ) { MyTitle=MyTitle." (".sChiDescr." ≤ ".sprintf("%g",chi_max).")" }
set title MyTitle
set ylabel 'a E_{eff}'
#set xlabel 'k^2 (sorted by '.sChiDescr.", ".num." best fits)"
set xlabel 'n^2 (sorted by '.sChiDescr.")"
set xtics 1
set key top left reverse Left
set pointsize 0.5

#plot for [i=0:MaxPSq] filename(i) index "[tf=".tf."]" every ::4::4 using (i):(\$3*\$3):(column(-2)) linecolor variable title "k^2=".i, \
    disp(x) lc rgb "gray20" lw 0.25  dashtype "-" title "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2"

#set label "Dispersion: (a E_0)^2 + (2 {/Times:Italic π} a / L)^2 k^2" at graph 1, graph 0 font "Arial,12" front textcolor "grey40" offset character -1.5, character 1.5 right
#set label "E_0=0.99656(95), ArXiv:1812.08791" at graph 1, graph 0 font "Arial,8" front textcolor "grey40" offset character -1, character 1 right
if( RefVal == -777 ) { RefText=MyField.'(0)='.gprintf("%g",E0).', ti='.DispTI.', tf='.DispTFLabel }
set label 1 RefText at graph 1, graph 0 font "Arial,8" front textcolor "grey40" noenhanced \
    offset character -1, character 1 right
set xrange [-0.05:MaxPSq+.99]
set samples 1000
MyScale=1./num < 0.08 ? 0.08 : 1./num
count = 0
last_idx=-1
plot \
  disp(x) lc rgb "gray20" lw 0.25 dashtype "-" title "Continuum disp", \
  LatDisp(x) lc rgb "#FF1493" lw 0.5 title "Lattice disp", \
  for [i=0:MaxPSq] filename(i) using \
    (@MyOffset):(@Condition @YField):(@Condition @YFieldLow):(@Condition (count=count+1, @YFieldHigh)):(@MySeries) \
    with yerrorbars notitle lc variable pt 6 pointsize 0.5, \
  for [idx=TIMin:TIMax] filename(0) index 0 using 1:(1/0):(1/0):(1/0) \
    with yerrorbars title SeriesName(idx) lc idx - TIMin + 1 pt 6, \
  for [i=0:MaxPSq] filename(i) using \
    (@MyOffset):(@Condition (count=count+1, @YFieldHigh)):"tfLabel" @WithLabels 0.5 notitle
}

EOFMark

fi
