#!/usr/bin/env bash

#set -x
set -e

###################################################

# User input

###################################################

Do=${Do-linear cubic shrink}
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
# at the end, also show the discretisation estimate

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
eval "set term tikz createstyle standalone ".PDFSize." font '\\\\small' header '\\\\usepackage{physics}'"

set xlabel '\$E_K$ / GeV'
set ylabel '$\abs{\delta f_{'.FLabel.'}}$ / \%' #^{D_s \rightarrow K}
set xrange[${xrange:-*:*}]

FieldFile(Field,Seq,FF)=CompareDir.'/'.Field.'_'.sprintf("%c",Seq+64).'_'.FF.'.txt'

array PlotTitles[22]
PlotTitles[1]='omit $\delta^{\left(\textrm{FV}\right)}$'
PlotTitles[2]='\$\Delta_0$ 313 MeV, $\Delta_+$ -5 MeV'
PlotTitles[3]='\$\Delta_+$ -5 MeV'
PlotTitles[4]='\$\Delta_0$ 313 MeV'
PlotTitles[5]='alternate \$Z_{\textrm{V,hh}}$'
PlotTitles[6]='alternate C1 fits'
PlotTitles[7]='omit \$n^2_{\textrm{max}}$ C,M'
PlotTitles[8]='omit C1'
PlotTitles[9]='omit C2'
PlotTitles[10]='omit M1'
PlotTitles[11]='omit M2'
PlotTitles[12]='omit M3'
PlotTitles[13]='omit M1M2M3'
PlotTitles[14]='\$Z_V\$ excited'
PlotTitles[15]='contin disp'

if( DoWhich eq 'cubic' ) {
  NumRows=22
# These are destructive tests - not part of my error budget
PlotTitles[16]='omit F1M'
PlotTitles[17]='omit C1C2 and $\left(a \Lambda\right)^2$'
PlotTitles[18]='omit $\delta^{\left(\textrm{chiral}\right)}$'
PlotTitles[19]=PlotTitles[1].' and $\delta^{\left(\textrm{chiral}\right)}$'
PlotTitles[20]='omit $\left(a \Lambda\right)^2$'
PlotTitles[21]='omit $\flatfrac{\Delta M_\pi^2}{\Lambda^2}$'
PlotTitles[22]='omit $\left(\flatfrac{E_L}{\Lambda}\right)^3$'
} else {
  NumRows=15
  #PlotTitles[16]='\$f_+ \left(\flatfrac{E_L}{\Lambda}\right)^3$ terms'
}


set linetype 5 lc rgb 0xFA60E4 #pink
set linetype 6 lc 'blue'
#set linetype 2 lc rgb 'dark-violet' # 0x9400D3
#set linetype 3 lc rgb 0xE36C09 #orange
set linetype 2 lc rgb 0x00B050 #green
#set linetype 5 lc rgb 0x0070C0 #blue
#set linetype 6 lc rgb 0x003860 #dark navy
#set linetype 7 lc rgb 0xFFC000 #yellow
#set linetype 8 lc rgb 0xC00000 #red
#set linetype 9 lc rgb 0xC00000 #black
NumColours=8

GetColour(i)=i<=7 ? i : (i-8)%NumColours+1
GetDashType(i)=i<=7 ? 'solid' : i<=13 ? '1' : '2'
set dashtype 1 (5, 5)
set dashtype 2 (1, 4)

PlotErrorBar="abs((column('y_high')-column('y_low'))/(2*column('y')))*100"
PlotPrefix="'$Ref/$Prefix'.FF.'$Suffix' using "
PlotPrefix=PlotPrefix."(stringcolumn('field') eq Field ? column('x')*XScale : NaN):("
PlotPrefixB="):("
PlotPrefixC=") with filledcurves title 'Stat error' fc 'gray05' fs transparent solid 0.5"

OutName=Save.'_'.FF

do for [Loop=1:2] {

ExtraName=''
if( Loop == 1 ) {
  ExtraName='_all'
}
if( Loop == 2 ) {
  NumRows=13
}

set key top left Left reverse maxrows (NumRows+3)/2 # NumRows + Error row and round up

AbsFunc='abs'
do for [AbsLoop=1:2] {

PlotPrefixA='0'
if( AbsLoop == 2 ) { PlotPrefixA='-'.PlotErrorBar }

Cmd=''
do for [i=1:NumRows] {
  if( PlotTitles[i] ne '' ) {
  Cmd=Cmd.', "'.FieldFile(Field,i,FF).'"'
  Cmd=Cmd." using (column('x')*XScale):(".AbsFunc."(column('Sys'))*100)"
  Cmd=Cmd." with lines title '".sprintf("%c",i+96).') '.PlotTitles[i]."'"
  Cmd=Cmd.' linewidth 1.5'
  Cmd=Cmd.' linetype '.GetColour(i)
  Cmd=Cmd.' dashtype '.GetDashType(i)
  #if( i > 16 ) { Cmd=Cmd.' dashtype 4' } else {
  #  if( i > 8 ) { Cmd=Cmd.' dashtype 2' } }
  }
}

Cmd=Cmd.', '.PlotPrefix.'column("DiscRelEr")*100) with lines title "'.sprintf("%c",NumRows+97).') discretisation"'
Cmd=Cmd.' linewidth 1.5'
Cmd=Cmd.' linetype '.GetColour(NumRows+1)
Cmd=Cmd.' dashtype '.GetDashType(NumRows+1)

set output PlotDir.'/'.OutName.AbsFunc.ExtraName.'.tex'

eval "plot ".PlotPrefix.PlotPrefixA.PlotPrefixB.PlotErrorBar.PlotPrefixC.Cmd
set output
AbsFunc=''
}
}
EOFMark

(
cp gnuplot*.* $PlotDir
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
  local LabelFirst="${Label:0:1}"
  if [ "${Label,,}" == "ref" ] || [[ "${LabelFirst,}" < "n" ]]
  then
    echo "$Out"' \\' >> "$SummaryFile"
  fi
  echo "$Out"' \\' >> "$SummaryAll"
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
    echo "${DefPrefix}FZQZVal{${QSqZero%(*}}"
    echo "${DefPrefix}FZQMaxVal{${QSqMax%(*}}"
    echo "${DefPrefix}FPQZVal{${QSq[0]%(*}}"
    echo "${DefPrefix}FPQMaxVal{${QSq[1]%(*}}"
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
      echo "x Ref Stat Sys Total StatSq SysSq TotalSq"
      # Join gets: x y Ref_low Ref Ref_high
      join -j 2 -o 1.2,1.6,2.5,2.6,2.7 \
        <(gawk -e "$GCmd" $File) <(gawk -e "$GCmd" $FileRef) \
      | awk -f <(cat - <<-'ENDGAWK'
		{ Stat=($5-$3)/(2*$4) }
		{ StatSq=Stat*Stat }
		{ Sys=($2-$4) / $4 }
  		{ SysSq=Sys*Sys }
		{ TotalSq=StatSq+SysSq }
		{ Total=sqrt(TotalSq) }
		{ print $1, $4, Stat, Sys, Total, StatSq, SysSq, TotalSq }
ENDGAWK
      )
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
  RefDir=${1:+$1/}
  RefSeries=$2
  RefSuffix=$3
  RefBase=$RefSeries$RefSuffix # This is the reference fit base name
  Ref=$RefDir$RefBase # This is the path to the reference fit
}

function MakeCommon()
{
  MakeOne b $RefDir/${RefBase}_PoleSV # Move scalar and vector poles
  MakeOne c $RefDir/${RefBase}_PoleV # Move vector pole
  MakeOne d $RefDir/${RefBase}_PoleS # Move scalar pole
  MakeOne e $RefDir/AltZV$RefSuffix # Alternate Z_{V, h-h}
  MakeOne f $RefDir/${RefBase}_AltC1 # Alternate C1
  MakeOne g Omit/${RefBase}_NMaxCM # Omit n^2_max from C & M ensembles
  MakeOne h Omit/${RefBase}-EnsC1 # Omit C1
  MakeOne i Omit/${RefBase}-EnsC2 # Omit C2
  MakeOne j Omit/${RefBase}-EnsM1 # Omit M1
  MakeOne k Omit/${RefBase}-EnsM2 # Omit M2
  MakeOne l Omit/${RefBase}-EnsM3 # Omit M3
  MakeOne m Omit/${RefBase}-EnsM # Omit M1, M2, M3
  MakeOne n $RefDir/Jan24$RefSuffix # Alternate Z_V multiple wall-separations and excited-states
  MakeOne o $RefDir/renormC$RefSuffix # Continuum dispersion relation
}

function MakeOriginal()
{
  MakeStats Ref "$Ref" # Reference value
  MakeOne a Omit/${RefSeries}E3-CVZ34 # Omit FV
  MakeCommon
  # These are destructive tests - not part of my error budget
  MakeOne p Omit/${RefBase}-EnsF1M # Omit F1M
  MakeOne q Omit/${RefSeries}E3D0-CZ34-EnsC # Omit C1, C2
  MakeOne r Omit/${RefSeries}E3-CXZ34 # Omit Chiral
  MakeOne s Omit/${RefSeries}E2-CVXZ3 # Omit FV and Chiral
  MakeOne t Omit/${RefSeries}E3D0-CZ34 # Omit a Lambda
  MakeOne u Omit/${RefSeries}E3-C1Z34 # Omit Delta Mpi / Lambda
  MakeOne v Omit/${RefSeries}E2-CZ3 # Omit (E/Lambda)^3
}

function MakeShrink()
{
  MakeStats Ref "$Ref" # Reference value
  MakeOne a Omit/${RefBase}-CV # Omit FV
  MakeCommon
  MakeOne p $RefDir/${RefSeries}E3-CZ34 # Original reference fit with cubic energy in f_+
}

###################################################

# Perform all the continuum fit variation comparisons

###################################################

function MakeAllHeader()
{
  if ((Narrow))
  then
    cat <<-"EOFMARK"
    \begin{tabular}{|c|r|c|l|l|l|l|l|}
    \hline
    & $t^2$ & $\nu$ & $t^2_\nu$ & p-H & $f_0 \lr{0}$ & $f_0 \lr{q^2_{\textrm{max}}}$ & $f_+ \lr{q^2_{\textrm{max}}}$ \\
    \hline
EOFMARK
  else
    cat <<-"EOFMARK"
    \begin{tabular}{|c|r|c|l|l|l|l|l|l|l|}
    \hline
    & $t^2$ & $\nu$ & $t^2_\nu$ & p-H & p-$\chi^2$ & $f_0 \lr{0}$ & $f_0 \lr{q^2_{\textrm{max}}}$ & $f_+ \lr{0}$ & $f_+ \lr{q^2_{\textrm{max}}}$ \\
    \hline
EOFMARK
  fi
}

function MakeAllFooter()
{
  cat <<-"EOFMARK"
  \hline
  \end{tabular}
EOFMARK
}

function MakeAll()
{
  mkdir -p "$CompareDir"
  local BaseName="$CompareDir/${Field}_Summary"
  local SummaryFile="$BaseName.txt"
  local SummaryAll="${BaseName}_all.txt"
  local SummaryTex="$BaseName.tex"
  echo "% Comparison: $Do" > "$SummaryTex"
  echo "% Comparison: $Do" > "$SummaryFile"
  echo "% Comparison: $Do" > "$SummaryAll"
  MakeAllHeader >> "$SummaryFile"
  MakeAllHeader >> "$SummaryAll"

  case ${Do,,} in
    cubic)             MakeOriginal;;
    shrink | linear)   MakeShrink;;
  esac

  MakeAllFooter >> "$SummaryFile"
  MakeAllFooter >> "$SummaryAll"
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
eval "set term tikz createstyle standalone ".PDFSize." font '\\\\small' header '\\\\usepackage{physics}'"

set xlabel '\$E_K$ / GeV'
set ylabel '$\abs{\delta f_{'.FLabel.'}}$ / \%' #^{D_s \rightarrow K}
set xrange[${xrange:-*:*}]
set arrow 1 from first 1.046578, graph 0 to first 1.046578, graph 1 \
  nohead front lc black lw 1 dashtype (3,6)
set arrow 2 from first 0.495644, graph 0 to first 0.495644, graph 1 \
  nohead front lc black lw 1 dashtype (10,5)

set key top center Left reverse maxrows 1 # outside
set output OutFile

MyX="(column('x')*XScale)"

plot InFile using @MyX:(0):(column('Total')*100) \
      with filledcurves title 'Total' fc 'gray05' fs transparent solid 0.5, \
  '' using @MyX:(column('Stat')*100) with lines lc "blue" lw 2 title 'Stat', \
  '' using @MyX:(abs(column('Sys'))*100) with lines lc "red" lw 2 title 'Sys'
#  , \
#  '' using @MyX:(column('fit')*100) with lines lc rgb 0xFA60E4 lw 2 dt (10,5) title 'Fit', \
#  '' using @MyX:(column('disc')*100) with lines lc rgb 0xE36C09 lw 2 dt (3,6) title 'Disc'
set output
EOFMark

(
  cp gnuplot*.* $PlotDir
  cd $PlotDir
  pdflatex "$FileBase"
)
}

###################################################

# Parameters: list of labels a, b, c, etc
# where each file has: x y Ref_low Ref Ref_high
# make the maximum of the y-column as statistical error
# NB: The ref values are the same in every file
# Optional:
#   Quad: Set to anything to add systematic errors in quadrature rather than maximum

###################################################

function MakeMax()
{
  local Label=($*)
  local i FileA FileB FileBase FileData ff IsSingle
  for ff in f0 fplus
  do
    FileA="$CompareDir/${Field}_${Label[0]}_${ff}.txt"
    FileBase="${Field}_Max_${ff}"
    FileData="$CompareDir/$FileBase.txt"
    IsSingle=1
    (( ${#Label[@]} > 1 )) && IsSingle=0
    for (( i=1-IsSingle; i < ${#Label[@]}; ++i ))
    do
      if (( i > 1 )); then
        FileA="$FileData.tmp"
        mv "$FileData" "$FileA"
      fi
      FileB="$CompareDir/${Field}_${Label[i]}_${ff}.txt"
      join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2,2.3,2.4,2.5,2.6,2.7,2.8 "$FileA" "$FileB" \
      | awk -v Quad=$((0${Quad+1})) -v IsSingle=$IsSingle -f <(cat - <<-'ENDGAWK'
		# Processing records of the form: x Ref Stat Sys Total StatSq SysSq TotalSq
		$1=="x" { print; next }
		{ if (IsSingle) { Sys=$4; SysSq=$7 }
            else { if (Quad) { SysSq=$7 + $14; Sys=sqrt(SysSq) }
			else { if ($14>$7) { Sys=$11; SysSq=$14 } else { Sys=$4; SysSq=$7 } } } }
		{ TotalSq=$6+SysSq }
		{ Total=sqrt(TotalSq) }
		{ print $1, $2, $3, Sys, Total, $6, SysSq, TotalSq }
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
      cubic)            SetRef Simul renorm E3-CZ34;;
      shrink | linear)  SetRef Simul renorm;;
    esac
    export xrange='0.48:*'
    for Field in EL
    do
      if [ -v Make ]; then MakeAll; fi
      # MakeMax a b c d e f g # Show maximum
      Quad= MakeMax a b e f # Add systematic errors in quadrature
      save=Var ff=f0 DoPlot
      save=Var ff=fplus DoPlot
    done
  )
  fi
done
