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

Exp2pt=3
RangeSnk=9_22
RangeSrc=8_23

Ensemble=${Ensemble:-F1M}
InBase=${InBase:-/Volumes/QCD/tursa/semilep/data}
PW=${PW:-g5P_g5P}
Seed=${Seed:-1835672416}
Analyse=${Analyse:-analyse}
Corr=${Corr:-corr}
Ratio=${Ratio:-ratio.manual}
Fit=${Fit:-fit}
MELFit=${MELFit:-MELFit}
LongPaths=0${Long+1}

FitWhat=${FitWhat:-quark}

Plot=${Plot:-Plot}
PlotDataFrom=${PlotDataFrom:-4} # how far in from source and sink to plot data points

############################################################

# Derived

############################################################

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

Fit2ptDir=${Fit}
[ "$Exp2pt" != 2 ] && Fit2ptDir=${Fit2ptDir}$Exp2pt
Fit2ptDir=$Fit2ptDir/2ptp2

SuffixModel=g5P_g5W.model.$Seed.h5
FitMesonSnk=$Fit2ptDir/$MesonSnk/${MesonSnk}_p2_${pSnk}.corr_${RangeSnk}_${RangeSnk}.${SuffixModel}
FitMesonSrc=$Fit2ptDir/$MesonSrc/${MesonSrc}_p2_${pSrc}.corr_${RangeSrc}_${RangeSrc}.${SuffixModel}

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

OutPart1=${qSnk}_${qSrc}_${Gamma}_p2_${pSrc}_ps2_${pSnk}
OutPart2=dt_$DeltaTAll
OutPart3=E${Exp2pt}_${RangeSnk}_${RangeSrc}
OutSubDir=$MELFit/$SpecDir
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

AllFitTimes=corr
for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  FitList="$FitList${MySep}${PrefixFit}${DeltaT[i]}${SuffixFit}.h5,t=${TI[i]}:${TF[i]}"
  DataList="$DataList${MySep}$DataDir/${PrefixFit}${DeltaT[i]}${SuffixFit}.txt"
  Title="$Title${MySep}ΔT=${DeltaT[i]}"
  LabelTI="$LabelTI${MySep}${PlotDataFrom}"
  LabelTF="$LabelTF${MySep}$(( ${DeltaT[i]} - PlotDataFrom ))"
  AllFitTimes=${AllFitTimes}_${TI[i]}_${TF[i]}
  MySep=" "
done

BuildModelBase=$OutSubDir/${FitWhat}_${OutLongName}

MultiFit="MultiFit -e 2 --Hotelling 0 --debug-signals --mindp 1 --overwrite"
Cmd="$MultiFit -i $DataDir/ -o $BuildModelBase $FitMesonSnk $FitMesonSrc $FitList"
BuildModelBase=$BuildModelBase.${AllFitTimes}
BuildModel=$BuildModelBase.g5P_g5W.model
echo "$Cmd"
echo "$Cmd"  > $BuildModel.$Seed.log
      $Cmd  >> $BuildModel.$Seed.log

# Plot it

mkdir -p $Plot/$OutSubDir
Cmd="files='$DataList' ti='$LabelTI' tf='$LabelTF' save=$Plot/${BuildModelBase}"
[ -v yrangeMEL ] && Cmd="$Cmd yrange='$yrangeMEL'"
Cmd="$Cmd title='$Title' field=$Field plottd.sh ${BuildModel}_td.$Seed.txt"
#echo "$Cmd"
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
PlotModelBase=$PlotModelBase.${AllFitTimes}
PlotModel=$PlotModelBase.g5P_g5W.model
echo "$Cmd"
echo "$Cmd"  > $PlotModel.$Seed.log
      $Cmd  >> $PlotModel.$Seed.log

# Plot it

mkdir -p $Plot/$OutSubDir
Cmd="files='$DataList' ti='$LabelTI' tf='$LabelTF' save=$Plot/${PlotModelBase}"
[ -v yrangeR3 ] && Cmd="$Cmd yrange='$yrangeR3'"
Cmd="$Cmd title='$Title' field=$PlotField plottd.sh ${PlotModel}_td.$Seed.txt"
#echo "$Cmd"
eval  $Cmd

