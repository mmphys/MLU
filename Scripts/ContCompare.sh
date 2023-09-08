#!/usr/bin/env bash

#set -x
set -e

###################################################

# User input

###################################################

Which=${Which:-all} # all / some
Base=Cont
CompareDir=Var
PlotDir=Plot
Ref=Simul/renorm-CZ34 # This is the reference fit

Prefix=F3_K_Ds.corr_
SuffixShort=.g5P_g5W.model
Suffix=${SuffixShort}_fit.txt
HName=${Prefix}f0_fplus.g5P_g5W.model

###################################################

# Plot all the continuum fit variation comparisons

###################################################

function DoPlot()
{
mkdir -p $PlotDir
gnuplot <<-EOFMark

PDFSize="${size:-6in,3in}"
Field="$Field"
PlotDir="$PlotDir"
Save="$save"
FF="$ff"
CompareDir="$CompareDir"

XScale=1e-9
if( Field eq "qSq" ) { XScale=XScale*XScale }
if( FF eq 'f0' ) { FLabel='0' } else { FLabel='+' }

if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
eval "set term tikz standalone ".PDFSize." font '\\\\small' header '\\\\usepackage{physics}'"

set xlabel '\$E_K$ / GeV'
set ylabel '$\delta f_{'.FLabel.'}^{D_s \rightarrow K}$ / \%'

FieldFile(Field,Seq,FF)=CompareDir.'/'.Field.'_'.sprintf("%c",Seq+64).'_'.FF.'.txt'

array PlotTitles[19]
PlotTitles[1]='omit $\delta^{\left(\textrm{FV}\right)}$'
PlotTitles[2]='omit C1'
PlotTitles[3]='omit C2'
PlotTitles[4]='omit M1'
PlotTitles[5]='omit M2'
PlotTitles[6]='omit M3'
PlotTitles[7]='omit M1M2M3'
PlotTitles[8]='omit \$n^2_{\textrm{max}}$ C1C2'
PlotTitles[9]='\$\Delta_0$ 300 MeV, $\Delta_+$ -25 MeV'
PlotTitles[10]='\$\Delta_+$ -100 MeV, omit $\left(\flatfrac{E_L}{\Lambda}\right)^3$'
PlotTitles[11]='\$\Delta_0$ 250 MeV'
PlotTitles[12]='Alternate C1'

# These are destructive tests - not part of my error budget
PlotTitles[13]='omit F1M'
PlotTitles[14]='omit C1C2 and $\left(a \Lambda\right)^2$'
PlotTitles[15]='omit $\delta^{\left(\textrm{chiral}\right)}$'
PlotTitles[16]=PlotTitles[1].' and $\delta^{\left(\textrm{chiral}\right)}$'
PlotTitles[17]='omit $\left(a \Lambda\right)^2$'
PlotTitles[18]='omit $\flatfrac{\Delta M_\pi^2}{\Lambda^2}$'
PlotTitles[19]='omit $\left(\flatfrac{E_L}{\Lambda}\right)^3$'

PlotPrefix="'$Ref/$Prefix'.FF.'$Suffix' using "
PlotPrefix=PlotPrefix."(stringcolumn('field') eq Field ? column('x')*XScale : NaN)"
PlotPrefix=PlotPrefix.":(0):(abs((column('y_high')-column('y_low'))/(2*column('y')))*100)"
PlotPrefix=PlotPrefix."with filledcurves title 'Error' fc 'gray05' fs transparent solid 0.5"

do for [Loop=1:2] {

if( Loop == 2 ) {
  do for [i=13:|PlotTitles|] { PlotTitles[i]='' }
  set key top left Left reverse maxrows 7
} else {
  set key center left Left reverse maxrows 10
}

Cmd=''
do for [i=1:|PlotTitles|] {
  if( PlotTitles[i] ne '' ) {
  Cmd=Cmd.', "'.FieldFile(Field,i,FF).'"'
  Cmd=Cmd." using (column('x')*XScale):(abs(column('y')/column('Ref')-1)*100)"
  Cmd=Cmd." with lines title '".sprintf("%c",i+96).') '.PlotTitles[i]."'"
  Cmd=Cmd.' linestyle '.i
  if( i > 16 ) { Cmd=Cmd.' dashtype 4' } else {
    if( i > 8 ) { Cmd=Cmd.' dashtype 2' } }
  }
}

OutName=Save.'_'.FF
if( Loop == 2 ) { OutName=OutName.'_zoom' }
set output PlotDir.'/'.OutName.'.tex'

eval "plot ".PlotPrefix.Cmd
set output
}
EOFMark

(
cd $PlotDir
pdflatex "${save}_${ff}"
pdflatex "${save}_${ff}_zoom"
)
}

###################################################

# Make statistics for a single fit

###################################################

function MakeStats()
{
  local Label="$1"
  local Dir=$2
  local File ff QSqZero QSqMax DoF ChiSq
  local -a ColumnValues QSq
  for ff in f0 fplus
  do
    File=$Dir/$Prefix$ff$Suffix
    # Now get stats
    QSq=($(gawk -e '/# qSq 0 / {print $4; next}; /# qSq / {print $4; exit}' $File))
    if [ "$ff" == "f0" ]; then
      QSqZero="${QSq[0]}"
      QSqMax="${QSq[1]}"
    fi
  done
  ColumnValues=($(GetColumn --exact ChiSqPerDof,pValueH,pValue $Dir/$HName.h5))
  DoF=($(h5dump -a model/DoF $Dir/$HName.h5 | gawk -e '/:/ {print $2; exit}'))
  ChiSq=$(bc <<< "scale=3;${ColumnValues[0*8+4]}*$DoF/1") # /1 to make sure scale works
  echo "$Label & $ChiSq & $DoF & ${ColumnValues[0*8+4]} & ${ColumnValues[1*8+4]} & ${ColumnValues[2*8+4]}" \
       "& $QSqZero & $QSqMax & ${QSq[0]} & ${QSq[1]} \\\\" \
    >> "$SummaryFile"
}

###################################################

# Make a single comparison between one fit and the reference file

###################################################

function MakeOne()
{
  local Label="$1"
  local Dir=$2
  local GCmd="\$1==\"$Field\""
  local File FileRef ff
  for ff in f0 fplus
  do
    File=$Dir/$Prefix$ff$Suffix
    FileRef=$Ref/$Prefix$ff$Suffix
    {
      echo "x y Ref"
      join -j 2 -o 1.2,1.6,2.6 <(gawk -e "$GCmd" $File) <(gawk -e "$GCmd" $FileRef)
    } > "$CompareDir/${Field}_${Label}_${ff}.txt"
  done
  MakeStats "$Label" "$Dir"
  (
    shopt -s nullglob
    for File in $Dir/${Prefix}*${SuffixShort}_pcorrel.txt; do
      echo "PlotMatrix.sh $File"
      PlotMatrix.sh "$File"
    done
  )
}

###################################################

# Perform all the continuum fit variation comparisons

###################################################

function MakeAll()
{
  mkdir -p "$CompareDir"
  local SummaryFile="$CompareDir/${Field}_Summary.txt"
  cat > "$SummaryFile" <<-"EOFMARK"
	\begin{tabular}{|c|r|c|l|l|l|l|l|l|l|}
	\hline
	& $\chi^2$ & DoF & $\chi^2$/dof & p-H & p-$\chi^2$ & $f_0 \lr{0}$ & $f_0 \lr{q^2_{\textrm{max}}}$ & $f_+ \lr{0}$ & $f_+ \lr{q^2_{\textrm{max}}}$ \\
	\hline
EOFMARK

  MakeStats Ref "$Ref" # Reference value
  MakeOne a Omit/renorm-CVZ34 # Omit FV
  MakeOne b C1/renorm-CZ34 # Omit C1
  MakeOne c C2/renorm-CZ34 # Omit C2
  MakeOne d M1/renorm-CZ34 # Omit M1
  MakeOne e M2/renorm-CZ34 # Omit M2
  MakeOne f M3/renorm-CZ34 # Omit M3
  MakeOne g CF/renorm-CZ34 # Omit M1, M2, M3
  MakeOne h Simul/renorm-CZ34_some # Omit n^2_max from C1 C2
  MakeOne i Pole300-25/renorm-CZ34 # Delta_+=-25MeV, Delta_0=300MeV
  MakeOne j PoleV-100/renormE2-CZ3 # Delta_+=-100MeV, Omit (E/Lambda)^3
  MakeOne k PoleS250/renorm-CZ34 # Delta_0=250MeV
  MakeOne l AltC1/renorm-CZ34 # Alternate C1

  # These are destructive tests - not part of my error budget
  MakeOne m F1M/renorm-CZ34 # Omit F1M
  MakeOne n MF/renormD0-CZ34 # Omit C1, C2
  MakeOne o Omit/renorm-CXZ34 # Omit Chiral
  MakeOne p Omit/renormE2-CVXZ3 # Omit FV and Chiral
  MakeOne q Omit/renormD0-CZ34 # Omit a Lambda
  MakeOne r Omit/renorm-C1Z34 # Omit Delta Mpi / Lambda
  MakeOne s Omit/renormE2-CZ3 # Omit (E/Lambda)^3

  cat >> "$SummaryFile" <<-"EOFMARK"
	\hline
	\end{tabular}
EOFMARK
}

###################################################

# Main loop

###################################################

if ! [ -d "$Base" ]; then echo "$Base doesn't exist"; exit 1; fi
cd "$Base"
if ! [ -d "$Ref" ]; then echo "Reference $Base/$Ref doesn't exist"; exit 1; fi

PlotMatrix.sh $Ref/${HName}_pcorrel.txt

for Field in EL
do
  if [ -v Make ]; then MakeAll; fi
  save=Var ff=f0 DoPlot
  save=Var ff=fplus DoPlot
done
