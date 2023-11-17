#!/usr/bin/env bash

#set -x
set -e

###################################################

# User input

###################################################

Do=${Do-cubic shrink linear}
Cont=Cont
CompareDir=Var
PlotDir=Plot

Prefix=F3_K_Ds.corr_
SuffixShort=.g5P_g5W.model
Suffix=${SuffixShort}_fit.txt
HName=${Prefix}f0_fplus.g5P_g5W.model
Narrow=$((1-0${Narrow+1})) # Narrow is the default (wasn't originally)

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
DoWhich="$Do"

XScale=1e-9
if( Field eq "qSq" ) { XScale=XScale*XScale }
if( FF eq 'f0' ) { FLabel='0' } else { FLabel='+' }

if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
eval "set term tikz standalone ".PDFSize." font '\\\\small' header '\\\\usepackage{physics}'"

set xlabel '\$E_K$ / GeV'
set ylabel '$\delta f_{'.FLabel.'}$ / \%' #^{D_s \rightarrow K}
set xrange[${xrange:-*:*}]

FieldFile(Field,Seq,FF)=CompareDir.'/'.Field.'_'.sprintf("%c",Seq+64).'_'.FF.'.txt'

array PlotTitles[19]
PlotTitles[1]='omit $\delta^{\left(\textrm{FV}\right)}$'
PlotTitles[2]='omit C1'
PlotTitles[3]='omit C2'
PlotTitles[4]='omit M1'
PlotTitles[5]='omit M2'
PlotTitles[6]='omit M3'
PlotTitles[7]='omit M1M2M3'
PlotTitles[8]='omit \$n^2_{\textrm{max}}$ C,M'
PlotTitles[9]='\$\Delta_0$ 313 MeV, $\Delta_+$ -5 MeV'
#PlotTitles[10]='\$\Delta_+$ -100 MeV, omit $\left(\flatfrac{E_L}{\Lambda}\right)^3$'
PlotTitles[10]='\$\Delta_+$ -5 MeV'
PlotTitles[11]='\$\Delta_0$ 313 MeV'
PlotTitles[12]='Alternate C1 fits'

if( DoWhich eq 'cubic' ) {
  NumRows=19
# These are destructive tests - not part of my error budget
PlotTitles[13]='omit F1M'
PlotTitles[14]='omit C1C2 and $\left(a \Lambda\right)^2$'
PlotTitles[15]='omit $\delta^{\left(\textrm{chiral}\right)}$'
PlotTitles[16]=PlotTitles[1].' and $\delta^{\left(\textrm{chiral}\right)}$'
PlotTitles[17]='omit $\left(a \Lambda\right)^2$'
PlotTitles[18]='omit $\flatfrac{\Delta M_\pi^2}{\Lambda^2}$'
PlotTitles[19]='omit $\left(\flatfrac{E_L}{\Lambda}\right)^3$'
} else {
  NumRows=13
  PlotTitles[13]='\$f_+ \left(\flatfrac{E_L}{\Lambda}\right)^3$ terms'
}

PlotErrorBar="abs((column('y_high')-column('y_low'))/(2*column('y')))*100"
PlotPrefix="'$Ref/$Prefix'.FF.'$Suffix' using "
PlotPrefix=PlotPrefix."(stringcolumn('field') eq Field ? column('x')*XScale : NaN):("
PlotPrefixB="):("
PlotPrefixC=") with filledcurves title 'Error' fc 'gray05' fs transparent solid 0.5"

OutName=Save.'_'.FF

do for [Loop=1:2] {

ExtraName=''
if( Loop == 1 ) {
  ExtraName='_all'
}
if( Loop == 2 ) {
  NumRows=12
}

set key top left Left reverse maxrows (NumRows+2)/2 # NumRows + Error row and round up

AbsFunc='abs'
do for [AbsLoop=1:2] {

PlotPrefixA='0'
if( AbsLoop == 2 ) { PlotPrefixA='-'.PlotErrorBar }

Cmd=''
do for [i=1:NumRows] {
  if( PlotTitles[i] ne '' ) {
  Cmd=Cmd.', "'.FieldFile(Field,i,FF).'"'
  Cmd=Cmd." using (column('x')*XScale):(".AbsFunc."(column('y')/column('Ref')-1)*100)"
  Cmd=Cmd." with lines title '".sprintf("%c",i+96).') '.PlotTitles[i]."'"
  Cmd=Cmd.' linestyle '.i
  if( i > 16 ) { Cmd=Cmd.' dashtype 4' } else {
    if( i > 8 ) { Cmd=Cmd.' dashtype 2' } }
  }
}

set output PlotDir.'/'.OutName.AbsFunc.ExtraName.'.tex'

eval "plot ".PlotPrefix.PlotPrefixA.PlotPrefixB.PlotErrorBar.PlotPrefixC.Cmd
set output
AbsFunc=''
}
}
EOFMark

(
cd $PlotDir
for Zoom in '' _all; do
  for AbsFunc in '' abs; do
    pdflatex "${save}_${ff}${AbsFunc}${Zoom}"
  done
done
)
}

###################################################

# Make statistics for a single fit

###################################################

function MakeStats()
{
  local Label="$1"
  local Dir=$2
  local Out File ff QSqZero QSqMax DoF ChiSq
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
  Out="$Label & $ChiSq & $DoF & ${ColumnValues[0*8+4]} & ${ColumnValues[1*8+4]}"
  if ((!Narrow)); then Out+=" & ${ColumnValues[2*8+4]}"; fi
  Out+=" & $QSqZero & $QSqMax"
  if ((!Narrow)); then Out+=" & ${QSq[0]}"; fi
  Out+=" & ${QSq[1]}"
  echo "$Out"' \\' >> "$SummaryFile"
  local DefPrefix='\def\ContFit'"$Label"
  {
    echo "${DefPrefix}ChiSq{$ChiSq}"
    echo "${DefPrefix}DoF{$DoF}"
    echo "${DefPrefix}ChiDoF{${ColumnValues[0*8+1]%(*}}"
    echo "${DefPrefix}PValH{${ColumnValues[1*8+1]%(*}}"
    echo "${DefPrefix}PVal{$(bc <<< "scale=4;(${ColumnValues[2*8+4]}+0.00005)/1")}"
    echo "${DefPrefix}FZQZ{$QSqZero}"
    echo "${DefPrefix}FZQMax{$QSqMax}"
    echo "${DefPrefix}FPQZ{${QSq[0]}}"
    echo "${DefPrefix}FPQMax{${QSq[1]}}"
  } >> "$SummaryTex"
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
      echo "x y Ref_low Ref Ref_high"
      join -j 2 -o 1.2,1.6,2.5,2.6,2.7 <(gawk -e "$GCmd" $File) <(gawk -e "$GCmd" $FileRef)
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

# Original and modified comparisons

###################################################

function SetRef()
{
  RefBase=$2 # This is the reference fit base name
  Ref=${1:+$1/}$RefBase # This is the path to the reference fit
}

function MakeCommon()
{
  MakeOne b Omit/${RefBase}-EnsC1 # Omit C1
  MakeOne c Omit/${RefBase}-EnsC2 # Omit C2
  MakeOne d Omit/${RefBase}-EnsM1 # Omit M1
  MakeOne e Omit/${RefBase}-EnsM2 # Omit M2
  MakeOne f Omit/${RefBase}-EnsM3 # Omit M3
  MakeOne g Omit/${RefBase}-EnsM # Omit M1, M2, M3
  MakeOne h Omit/${RefBase}_NMaxCM # Omit n^2_max from C & M ensembles
  MakeOne i Simul/${RefBase}_PoleSV # Move scalar and vector poles
  MakeOne j Simul/${RefBase}_PoleV # Move vector pole
  MakeOne k Simul/${RefBase}_PoleS # Move scalar pole
  MakeOne l Simul/${RefBase}_AltC1 # Alternate C1
}

function MakeOriginal()
{
  MakeStats Ref "$Ref" # Reference value
  MakeOne a Omit/renormE3-CVZ34 # Omit FV
  MakeCommon
  # These are destructive tests - not part of my error budget
  MakeOne m Omit/${RefBase}-EnsF1M # Omit F1M
  MakeOne n Omit/renormE3D0-CZ34-EnsC # Omit C1, C2
  MakeOne o Omit/renormE3-CXZ34 # Omit Chiral
  MakeOne p Omit/renormE2-CVXZ3 # Omit FV and Chiral
  MakeOne q Omit/renormE3D0-CZ34 # Omit a Lambda
  MakeOne r Omit/renormE3-C1Z34 # Omit Delta Mpi / Lambda
  MakeOne s Omit/renormE2-CZ3 # Omit (E/Lambda)^3
}

function MakeShrink()
{
  MakeStats Ref "$Ref" # Reference value
  MakeOne a Omit/${RefBase}-CV # Omit FV
  MakeCommon
  MakeOne m Simul/renormE3-CZ34 # Original reference fit with cubic energy in f_+
}

###################################################

# Perform all the continuum fit variation comparisons

###################################################

function MakeAll()
{
  mkdir -p "$CompareDir"
  local BaseName="$CompareDir/${Field}_Summary"
  local SummaryFile="$BaseName.txt"
  local SummaryTex="$BaseName.tex"
  echo "% Comparison: $Do" > "$SummaryTex"
  echo "% Comparison: $Do" > "$SummaryFile"
  if ((Narrow))
  then
    cat >> "$SummaryFile" <<-"EOFMARK"
    \begin{tabular}{|c|r|c|l|l|l|l|l|}
    \hline
    & $t^2$ & $\nu$ & $t^2_\nu$ & p-H & $f_0 \lr{0}$ & $f_0 \lr{q^2_{\textrm{max}}}$ & $f_+ \lr{q^2_{\textrm{max}}}$ \\
    \hline
EOFMARK
  else
    cat >> "$SummaryFile" <<-"EOFMARK"
    \begin{tabular}{|c|r|c|l|l|l|l|l|l|l|}
    \hline
    & $t^2$ & $\nu$ & $t^2_\nu$ & p-H & p-$\chi^2$ & $f_0 \lr{0}$ & $f_0 \lr{q^2_{\textrm{max}}}$ & $f_+ \lr{0}$ & $f_+ \lr{q^2_{\textrm{max}}}$ \\
    \hline
EOFMARK
  fi

  case ${Do,,} in
    cubic)             MakeOriginal;;
    shrink | linear)   MakeShrink;;
  esac

  cat >> "$SummaryFile" <<-"EOFMARK"
	\hline
	\end{tabular}
EOFMARK
}

###################################################

# Plot all the continuum fit variation comparisons

###################################################

function PlotMax()
{
mkdir -p $PlotDir
gnuplot <<-EOFMark

PDFSize="${size:-6in,2in}"
Field="$Field"
FF="$ff"
InFile="$FileData"
OutFile="$PlotDir/$FileBase.tex"

XScale=1e-9
if( Field eq "qSq" ) { XScale=XScale*XScale }
if( FF eq 'f0' ) { FLabel='0' } else { FLabel='+' }

if( PDFSize ne "" ) { PDFSize='size '.PDFSize }
eval "set term tikz standalone ".PDFSize." font '\\\\small' header '\\\\usepackage{physics}'"

set xlabel '\$E_K$ / GeV'
set ylabel '$\delta f_{'.FLabel.'}$ / \%' #^{D_s \rightarrow K}
set xrange[${xrange:-*:*}]
set arrow 1 from first 1.046578, graph 0 to first 1.046578, graph 1 \
  nohead front lc rgb "gray40" lw 0.25 dashtype ".  "
set arrow 2 from first 0.495644, graph 0 to first 0.495644, graph 1 \
  nohead front lc rgb "gray40" lw 0.25 dashtype "-  "

set key top center Left reverse maxrows 1
set output OutFile

plot InFile using (column('x')*XScale):(0):(column('stat')*100) with filledcurves lc "skyblue" fs transparent solid 0.5 title 'Stat', \
  '' using (column('x')*XScale):(column('stat')*100):(column('total')*100) with filledcurves lc "red" fs transparent solid 0.5 title 'Sys'

set output
EOFMark

(
  cd $PlotDir
  pdflatex "$FileBase"
)
}

###################################################

# Parameters: list of labels a, b, c, etc
# where each file has: x y Ref_low Ref Ref_high
# make the maximum of the y-column as statistical error
# NB: The ref values are the same in every file

###################################################

function MakeMax()
{
  local Label=($*)
  local i FileA FileB FileBase FileData ff
  for ff in f0 fplus
  do
    FileA="$CompareDir/${Field}_${Label[0]}_${ff}.txt"
    FileBase="${Field}_Max_${ff}"
    FileData="$CompareDir/$FileBase.txt"
    i=0
    (( ${#Label[@]} > 1 )) && i=1
    for (( ; i < ${#Label[@]}; ++i ))
    do
      if (( i > 1 )); then
        FileA="$FileData.tmp"
        mv "$FileData" "$FileA"
      fi
      FileB="$CompareDir/${Field}_${Label[i]}_${ff}.txt"
      join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 "$FileA" "$FileB" \
      | awk -f <(cat - <<-'ENDGAWK'
		#{ print }
		$1=="x" { print "x y Ref_low Ref Ref_high stat sys total"; next }
		{ Stat=($5-$3)/(2*$4) }
  		{ RelA=($2>$4 ? $2-$4 : $4-$2) / $4 }
		{ RelB=($6>$4 ? $6-$4 : $4-$6) / $4 }
		{ Max=RelA > RelB ? $2 : $6 }
		{ Sys=RelA > RelB ? RelA : RelB }
		{ print $1, Max, $3, $4, $5, Stat, Sys, Stat+Sys }
ENDGAWK
      ) > "$FileData"
      (( i > 1 )) && rm "$FileA"
    done
    PlotMax
  done
}

###################################################

# Main loop

###################################################

for Do in ${Do,,}
do
  if ! [ -d "$Cont/$Do" ]; then
    echo "Reference $Cont/$Do doesn't exist"
  else
  (
    cd "$Cont/$Do"
    case "$Do" in
      cubic)            SetRef Simul renormE3-CZ34;;
      shrink | linear)  SetRef Simul renorm;;
    esac
    export xrange='0.48:*'
    for Field in EL
    do
      if [ -v Make ]; then MakeAll; fi
      MakeMax a b c d e f g h i j k l
      save=Var ff=f0 DoPlot
      save=Var ff=fplus DoPlot
    done
  )
  fi
done
