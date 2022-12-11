#!/usr/bin/env bash

# Perform a fit for matrix elements given a choice of 2pt fit ranges

############################################################

# Inputs

############################################################

#Ensemble=${Ensemble:-F1M}
. PlotCommon.sh

#set -x
set -e

qSrc=${qSrc:-h$Heavy}
qSnk=${qSnk:-l}
qSpec=${qSpec:-s}

pSnk=${pSnk:-0}
pSrc=0
Gamma=${Gamma:-gT}

# Ranges for 2pt fits
FitSnk=${FitSnk:-corr_10_26_7_26}
FitSrc=${FitSrc:-corr_6_29_5_29}

DeltaT=(${DeltaT:-24 28 32})
TI=(${TI:-13 14 14})
TF=(${TF:-16 20 24})

case "$FitWhat" in
  R3 | quark) : ;;
          "") FitWhat=R3;;
           *) echo "FitWhat=$FitWhat unrecognised"; exit 1;;
esac
NumExp=${NumExp:-3} # How many exponentials in the 3pt fit

PW=${PW:-g5P_g5P}
Seed=${Seed:-1835672416}
Ratio=${Ratio:-ratioE1ZV1}
Corr=${Corr:-corr}
MELFit=${MELFit:-MELFit}
#yrangeR3=${yrangeR3:-0.00068:0.000782}
#yrangeMEL=${yrangeMEL:-0.527:0.537}

Plot=${Plot:-Plot}
PlotDataFrom=${PlotDataFrom:-4} # how far in from source and sink to plot data points

############################################################

# Functions

############################################################

############################################################

# Derived

############################################################

DeltaTAll=${DeltaT[@]}
DeltaTAll=${DeltaTAll// /_}

GetMesonFile MesonSnk $qSnk $qSpec
GetMesonFile MesonSrc $qSrc $qSpec

SpecDir=3sm_${qSpec}p2

PrefixR3=$PlotData/$Ratio/$SpecDir/R3_${qSnk}_${qSrc}_${Gamma}_dt_
SuffixR3=_p2_${pSnk}_${PW}.fold.${Seed}
PrefixCorr=$PlotData/$Corr/$SpecDir/quark_${qSnk}_${qSrc}_${Gamma}_dt_
SuffixCorr=_p2_${pSrc}_ps2_${pSnk}_${PW}.fold.${Seed}

case "$FitWhat" in
  quark)
        FitName="C^{(3)}"
        PlotName=R_3
        PlotOptions=",raw"
        PrefixFit=$PrefixCorr
        SuffixFit=$SuffixCorr
        PrefixPlot=$PrefixR3
        SuffixPlot=$SuffixR3
        PlotWhat=R3
        Field=log
        PlotField=corr
        [ -v yrangeMEL ] && yrangeFit="$yrangeMEL"
        [ -v yrangeR3 ] && yrangePlot="$yrangeR3";;
  R3)
        FitName=R_3
        PlotName="C^{(3)}"
        FitOptions=",raw"
        PrefixFit=$PrefixR3
        SuffixFit=$SuffixR3
        PrefixPlot=$PrefixCorr
        SuffixPlot=$SuffixCorr
        PlotWhat=quark
        Field=corr
        PlotField=log
        [ -v yrangeR3 ] && yrangeFit="$yrangeR3"
        [ -v yrangeMEL ] && yrangePlot="$yrangeMEL";;
  *) echo "Error: FitWhat=$FitWhat"; exit 1;;
esac

FitType=corr
SuffixModel=g5P_g5W.model
Suffix2pt=$SuffixModel.$Seed.h5
InputMesonSnk=$Ensemble/$MELFit/2ptp2/$MesonSnk/${MesonSnk}_p2_${pSnk}.$FitSnk.${Suffix2pt}
InputMesonSrc=$Ensemble/$MELFit/2ptp2/$MesonSrc/${MesonSrc}_p2_${pSrc}.$FitSrc.${Suffix2pt}

OutPart1=${qSnk}_${qSrc}_${Gamma}_p2_$((pSrc > pSnk ? pSrc : pSnk))
OutPart2=dt_$DeltaTAll
OutPart3=$FitSnk.$FitSrc.E$NumExp
OutSubDir=$Ensemble/$MELFit/$SpecDir
OutLongName=${OutPart1}.${OutPart2}.${OutPart3}

############################################################

# Main loop

############################################################

# Fit the data

mkdir -p $OutSubDir

for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  FitList="$FitList${MySep}${PrefixFit}${DeltaT[i]}${SuffixFit}.h5,t=${TI[i]}:${TF[i]}${FitOptions}"
  DataList="$DataList${DataList:+ }${PrefixFit}${DeltaT[i]}${SuffixFit}.txt"
  Title="$Title${MySep}ΔT=${DeltaT[i]}"
  LabelTI="$LabelTI${MySep}${PlotDataFrom}"
  LabelTF="$LabelTF${MySep}$(( ${DeltaT[i]} - PlotDataFrom ))"
  FitType=${FitType}_${TI[i]}_${TF[i]}
  MySep=" "
done

BuildModelBase=$OutSubDir/${FitWhat}_${OutLongName}

MultiFit="MultiFit -e $NumExp --Hotelling 0 --mindp 1"
MultiFit="$MultiFit --overwrite"
MultiFit="$MultiFit --debug-signals"
Cmd="$MultiFit --summary 2 -o $BuildModelBase $InputMesonSnk $InputMesonSrc $FitList"
BuildModelBase=$BuildModelBase.${FitType}
BuildModel=$BuildModelBase.$SuffixModel
#echo "A: $Cmd"
echo "$Cmd"  > $BuildModel.$Seed.log
if  ! $Cmd &>> $BuildModel.$Seed.log
then
  LastError=${PIPESTATUS[0]}
  if [ "$LastError" = 3 ]; then
    echo "Warning: Not all parameters resolved"
  else
    echo "Warning $LastError: $Cmd"
  fi
fi

# Get the fit characteristics: energy difference, matrix element, test stat, ...
ColumnValues=$(GetColumn --exact ChiSqPerDof,pValueH --partial EDiff,MEL0,R3Raw ${BuildModel}.$Seed.h5)
if [ "$?" != 0 ]; then
  echo "Error: $ColumnValues"
  unset ColumnValues
  unset RefText
else
  #echo "OK: $ColumnValues"
  ColumnValues=($ColumnValues)
  EDiff="${ColumnValues[@]:16:8}"
  MEL="${ColumnValues[@]:24:8}"
  R3Raw="${ColumnValues[@]:32:8}"
  RefText="MEL0=${ColumnValues[25]} χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})"
fi

case "$FitWhat" in
  quark)
        [ -v EDiff ] && RefValFit="$EDiff"
        [ -v R3Raw ] && RefValPlot="$R3Raw";;
  R3)
        [ -v R3Raw ] && RefValFit="$R3Raw"
        [ -v EDiff ] && RefValPlot="$EDiff";;
esac

# Plot it

Cmd="files='$DataList' ti='$LabelTI' tf='$LabelTF' save=${BuildModelBase} $ModelMin"
[ -v yrangeFit ] && Cmd="$Cmd yrange='$yrangeFit'"
[ -v RefText ] && Cmd="$Cmd RefText='Fit $FitName: $RefText'"
[ -v RefValFit ] && Cmd="$Cmd RefVal='$RefValFit'"
Cmd="$Cmd title='$Title' field=$Field plottd.sh ${BuildModel}_td.$Seed.txt"
#echo "B: $Cmd"
echo "$Cmd"   > $BuildModelBase.log
eval  $Cmd  &>> $BuildModelBase.log

# Use the model to create the alternate

OutSubDirName=$OutSubDir/2Plot

unset MySep
unset FitList
unset DataList
for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  ThisFile=${PrefixPlot}${DeltaT[i]}${SuffixPlot}
  FitList="$FitList${MySep}${ThisFile}.h5,t=${TI[i]}:${TF[i]}${PlotOptions}"
  DataList="$DataList${MySep}${ThisFile}.txt"
  MySep=" "
done

mkdir -p $OutSubDirName
PlotModelBase=$OutSubDirName/${PlotWhat}_${OutLongName}

Cmd="$MultiFit -o $PlotModelBase ${BuildModel}.$Seed.h5 $FitList"
PlotModelBase=$PlotModelBase.${FitType}
PlotModel=$PlotModelBase.$SuffixModel
#echo "C: $Cmd"
echo "$Cmd"  > $PlotModel.$Seed.log
if !  $Cmd &>> $PlotModel.$Seed.log
then
  LastError=${PIPESTATUS[0]}
  echo "Warning $LastError: $Cmd"
fi

# Plot it

Cmd="files='$DataList' ti='$LabelTI' tf='$LabelTF' save=${PlotModelBase}"
[ -v yrangePlot ] && Cmd="$Cmd yrange='$yrangePlot'"
[ -v RefText ] && Cmd="$Cmd RefText='Plot $PlotName: $RefText'"
[ -v RefValPlot ] && Cmd="$Cmd RefVal='$RefValPlot'"
Cmd="$Cmd title='$Title' field=$PlotField plottd.sh ${PlotModel}_td.$Seed.txt"
#echo "D: $Cmd"
echo "$Cmd"      > $PlotModelBase.log
if ! eval $Cmd &>> $PlotModelBase.log
then
  LastError=${PIPESTATUS[0]}
  echo "Error $LastError: $Cmd"
fi
