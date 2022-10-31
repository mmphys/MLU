#!/usr/bin/env bash
#. common_utility.sh
. PlotCommon.sh

#set -x
set -e

############################################################

# Inputs

############################################################

qSrc=${qSrc:-h385}
qSnk=${qSnk:-l}
qSpec=${qSpec:-s}

pSnk=${pSnk:-0}
pSrc=0
Gamma=gT

DeltaT=(${DeltaT:-20 24 28 32})
TI=(${TI:-10 12 14 10})
TF=(${TF:-13 16 20 26})

# Ranges for 2pt fits
Exp2pt=${Exp2pt:-3}
TISnk=${TISnk:-9}
TFSnk=${TFSnk:-22}
TISrc=${TISrc:-8}
TFSrc=${TFSrc:-23}
Simul=0${Simul+1} # Set to anything to perform simultaneous 2- and 3-pt fit

Ensemble=${Ensemble:-F1M}
InBase=${InBase:-/Volumes/QCD/tursa/semilep/data}
PW=${PW:-g5P_g5P}
Seed=${Seed:-1835672416}
Analyse=${Analyse:-analyse}
Corr=${Corr:-corr}
Ratio=${Ratio:-ratio}
Fit=${Fit:-fit}
MELFit=${MELFit:-MELFit}
LongPaths=0${Long+1}

FitWhat=${FitWhat:-quark}

Plot=${Plot:-Plot}
PlotDataFrom=${PlotDataFrom:-4} # how far in from source and sink to plot data points

############################################################

# Derived

############################################################

# Ranges for 2pt fits
RangeSnk=${TISnk}_${TFSnk}
RangeSrc=${TISrc}_${TFSrc}

DeltaTAll=${DeltaT[@]}
DeltaTAll=${DeltaTAll// /_}

GetMesonFile MesonSnk $qSnk $qSpec
GetMesonFile MesonSrc $qSrc $qSpec

SpecDir=3sm_${qSpec}p2

DataDir=$InBase/$Ensemble/$Analyse

PrefixR3=R3_${qSnk}_${qSrc}_${Gamma}_dt_
SuffixR3=_p2_${pSnk}_${PW}.fold.${Seed}
PrefixCorr=quark_${qSnk}_${qSrc}_${Gamma}_dt_
SuffixCorr=_p2_${pSrc}_ps2_${pSnk}_${PW}.fold.${Seed}

case "$FitWhat" in
  quark)
        PrefixFit=$Corr/$SpecDir/$PrefixCorr
        SuffixFit=$SuffixCorr
        PrefixPlot=$Ratio/$SpecDir/$PrefixR3
        SuffixPlot=$SuffixR3
        PlotWhat=R3
        Field=log
        PlotField=corr;;
  R3)
        PrefixFit=$Ratio/$SpecDir/$PrefixR3
        SuffixFit=$SuffixR3
        PrefixPlot=$Corr/$SpecDir/$PrefixCorr
        SuffixPlot=$SuffixCorr
        PlotWhat=quark
        Field=corr
        PlotField=log;;
  *) echo "Error: FitWhat=$FitWhat"; exit 1;;
esac

FitType=corr
FitTypePlot=$FitType
if (( Simul ))
then
  FitType=${FitType}_${RangeSnk}_${RangeSrc}
  Input2ptDir=${Corr}/2ptp2
  SuffixModel=g5P.model
  Suffix2pt=g5P_g5P.fold.$Seed
  InputMesonSnk=$Input2ptDir/${MesonSnk}_p2_${pSnk}_${Suffix2pt}
  InputMesonSrc=$Input2ptDir/${MesonSrc}_p2_${pSrc}_${Suffix2pt}
  DataList="$DataDir/$InputMesonSnk.txt $DataDir/$InputMesonSrc.txt"
  InputMesonSnk=$InputMesonSnk.h5,t=${TISnk}:${TFSnk},e=${Exp2pt}
  InputMesonSrc=$InputMesonSrc.h5,t=${TISrc}:${TFSrc},e=${Exp2pt}
  OutPart3=S
  ModelMin="mmin=2"
else
  Input2ptDir=${Fit}
  [ "$Exp2pt" != 2 ] && Input2ptDir=${Input2ptDir}$Exp2pt
  Input2ptDir=$Input2ptDir/2ptp2
  SuffixModel=g5P_g5W.model
  Suffix2pt=$SuffixModel.$Seed.h5
  InputMesonSnk=$Input2ptDir/$MesonSnk/${MesonSnk}_p2_${pSnk}.corr_${RangeSnk}_${RangeSnk}.${Suffix2pt}
  InputMesonSrc=$Input2ptDir/$MesonSrc/${MesonSrc}_p2_${pSrc}.corr_${RangeSrc}_${RangeSrc}.${Suffix2pt}
  OutPart3=E
fi

OutPart1=${qSnk}_${qSrc}_${Gamma}_p2_${pSrc}_ps2_${pSnk}
OutPart2=dt_$DeltaTAll
OutPart3=${OutPart3}${Exp2pt}_${RangeSnk}_${RangeSrc}
OutSubDir=$MELFit/$Ensemble/$SpecDir
((LongPaths)) && OutSubDir=$OutSubDir/$OutPart1/$OutPart2/$OutPart3/${TI}_${TF}
OutLongName=${OutPart1}.${OutPart2}.${OutPart3}

############################################################

# Functions

############################################################

############################################################

# Main loop

############################################################

# Fit the data

mkdir -p $OutSubDir

for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  FitList="$FitList${MySep}${PrefixFit}${DeltaT[i]}${SuffixFit}.h5,t=${TI[i]}:${TF[i]}"
  DataList="$DataList${DataList:+ }$DataDir/${PrefixFit}${DeltaT[i]}${SuffixFit}.txt"
  Title="$Title${MySep}ΔT=${DeltaT[i]}"
  LabelTI="$LabelTI${MySep}${PlotDataFrom}"
  LabelTF="$LabelTF${MySep}$(( ${DeltaT[i]} - PlotDataFrom ))"
  FitType=${FitType}_${TI[i]}_${TF[i]}
  FitTypePlot=${FitTypePlot}_${TI[i]}_${TF[i]}
  MySep=" "
done

BuildModelBase=$OutSubDir/${FitWhat}_${OutLongName}

MultiFit="MultiFit -e 2 --Hotelling 0 --mindp 1"
MultiFit="$MultiFit --overwrite"
MultiFit="$MultiFit --debug-signals"
Cmd="$MultiFit --summary 2 -i $DataDir/ -o $BuildModelBase $InputMesonSnk $InputMesonSrc $FitList"
BuildModelBase=$BuildModelBase.${FitType}
BuildModel=$BuildModelBase.$SuffixModel
echo "A: $Cmd"
echo "$Cmd"  #> $BuildModel.$Seed.log
      $Cmd  >> $BuildModel.$Seed.log

# Get the energy difference
ColumnValues=$(GetColumn --exact ChiSqPerDof --partial EDiff,MEL0 ${BuildModel}.$Seed.h5)
if [ "$?" != 0 ]; then
  echo "Error: $ColumnValues"
  unset ColumnValues
else
  echo "OK: $ColumnValues"
  ColumnValues=($ColumnValues)
  ChiSqPerDof="χ²/dof=${ColumnValues[@]:3:1}"
  EDiff="${ColumnValues[@]:7:7}"
  MEL="${ColumnValues[@]:14:7}"
fi

# Plot it

mkdir -p $Plot/$OutSubDir
Cmd="files='$DataList' ti='$LabelTI' tf='$LabelTF' save=$Plot/${BuildModelBase} $ModelMin"
[ -v yrangeMEL ] && Cmd="$Cmd yrange='$yrangeMEL'"
[ -v EDiff ] && [ "$Field" = log ] && Cmd="$Cmd RefVal='$EDiff'"
[ -v MEL ] && [ "$Field" = corr ] && Cmd="$Cmd RefVal='$MEL'"
[ -v ChiSqPerDof ] && Cmd="$Cmd RefText='$ChiSqPerDof'"
Cmd="$Cmd title='$Title' field=$Field plottd.sh ${BuildModel}_td.$Seed.txt"
#echo "B: $Cmd"
eval  $Cmd

# Use the model to create the alternate

unset MySep
unset FitList
unset DataList
for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  ThisFile=$DataDir/${PrefixPlot}${DeltaT[i]}${SuffixPlot}
  FitList="$FitList${MySep}${ThisFile}.h5,t=${TI[i]}:${TF[i]}"
  DataList="$DataList${MySep}${ThisFile}.txt"
  MySep=" "
done

PlotModelBase=$OutSubDir/${PlotWhat}_${OutLongName}

Cmd="$MultiFit -o $PlotModelBase ${BuildModel}.$Seed.h5 $FitList"
PlotModelBase=$PlotModelBase.${FitTypePlot}
PlotModel=$PlotModelBase.$SuffixModel
#echo "C: $Cmd"
echo "$Cmd"  > $PlotModel.$Seed.log
      $Cmd  >> $PlotModel.$Seed.log

# Plot it

mkdir -p $Plot/$OutSubDir
Cmd="files='$DataList' ti='$LabelTI' tf='$LabelTF' save=$Plot/${PlotModelBase}"
[ -v yrangeR3 ] && Cmd="$Cmd yrange='$yrangeR3'"
[ -v EDiff ] && [ "$PlotField" = log ] && Cmd="$Cmd RefVal='$EDiff'"
[ -v MEL ] && [ "$PlotField" = corr ] && Cmd="$Cmd RefVal='$MEL'"
[ -v ChiSqPerDof ] && Cmd="$Cmd RefText='$ChiSqPerDof'"
Cmd="$Cmd title='$Title' field=$PlotField plottd.sh ${PlotModel}_td.$Seed.txt"
#echo "D: $Cmd"
eval  $Cmd

