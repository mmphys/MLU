#!/usr/bin/env bash

############################################################

# Functions for performing various versions of matrix element fits

############################################################

. PlotCommon.sh

#set -x
set -e

############################################################

# Perform two-stage fit

############################################################

# Parameters:
#   1: pSnk
function FitTwoStage()
(
  local pSnk=$1
  local MesonSnk; GetMesonFile MesonSnk $qSnk $qSpec
  local MesonSrc; GetMesonFile MesonSrc $qSrc $qSpec
  local FitSnk=${aMesonFit[$MesonSnk,$pSnk]}
  local FitSrc=${aMesonFit[$MesonSrc,0]}
  local FileOpSnk=${aMesonFileOp[$MesonSnk,$pSnk]}
  local FileOpSrc=${aMesonFileOp[$MesonSrc,0]}
  local FileMomSnk=${aMesonFileMom[$MesonSnk,$pSnk]}
  local FileMomSrc=${aMesonFileMom[$MesonSrc,0]}

############################################################

# Inputs

############################################################

qSrc=${qSrc:-h$Heavy}
qSnk=${qSnk:-l}
qSpec=${qSpec:-s}

pSnk=${pSnk:-0}
pSrc=0
Gamma=${Gamma:-gT}

# Ranges for 2pt fits
FitSnk=${FitSnk:-corr_10_26_7_26}
FitSrc=${FitSrc:-corr_6_29_5_29}
FileOpSnk=${FileOpSnk:-g5P_g5W}
FileOpSrc=${FileOpSrc:-g5P_g5W}
FileMomSnk=${FileMomSnk-_p2_$pSnk}
FileMomSrc=${FileMomSrc-_p2_$pSrc}

DeltaT=(${DeltaT:-24 28 32})
TI=(${TI:-13 14 14})
TF=(${TF:-16 20 24})

case "$FitWhat" in
  R3 | quark) : ;;
          "") FitWhat=R3;;
           *) echo "FitWhat=$FitWhat unrecognised"; exit 1;;
esac
NumExp=${NumExp:-3} # How many exponentials in the 3pt fit
FileSeries=${FileSeries:-E$NumExp}

PW=${PW:-g5P_g5P}
MLUSeed=${MLUSeed:-1835672416}
Ratio=${Ratio:-ratioE1ZV1}
Corr=${Corr:-corr}
MELFit=${MELFit:-MELFit}
yrangeR3=${yrangeR3:-${ayrangeR3[$Gamma,$pSnk]}}
yrangeMEL=${yrangeMEL:-${ayrangeMEL[$Gamma,$pSnk]}}
eval FitOptionsPerFile=($Options)
eval ThinningPerFile=($Thinning)

Plot=${Plot:-Plot}
PlotDataFrom=${PlotDataFrom:-4} # how far in from source and sink to plot data points

# UnCorr: set to anything to perform uncorrelated fit
# FitOptions: extra options for MultiFit command-line

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
SuffixR3=_p2_${pSnk}_${PW}.fold.${MLUSeed}
PrefixCorr=$PlotData/$Corr/$SpecDir/quark_${qSnk}_${qSrc}_${Gamma}_dt_
SuffixCorr=_p2_${pSrc}_ps2_${pSnk}_${PW}.fold.${DataSeed}

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
        FitOptionsModel=",C2e=2,C2Model=cosh,raw"
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

[ -v UnCorr ] && FitType=uncorr || FitType=corr
SuffixModel=g5P_g5W.model
MesonDir=$Ensemble/$MELFit/2ptp2
InputMesonSnk=$MesonDir/$MesonSnk/${MesonSnk}${FileMomSnk}.$FitSnk.${FileOpSnk}.model.$MLUSeed.h5
InputMesonSrc=$MesonDir/$MesonSrc/${MesonSrc}${FileMomSrc}.$FitSrc.${FileOpSrc}.model.$MLUSeed.h5

OutPart1=${qSnk}_${qSrc}_${Gamma}_p2_$((pSrc > pSnk ? pSrc : pSnk))
OutPart2=dt_$DeltaTAll
OutPart3=$FitSnk.$FitSrc.$FileSeries
OutSubDir=$Ensemble/$MELFit/$SpecDir
OutLongName=${OutPart1}.${OutPart2}.${OutPart3}

############################################################

# Main loop

############################################################

# Fit the data

mkdir -p $OutSubDir

for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  FitList="$FitList${MySep}${PrefixFit}${DeltaT[i]}${SuffixFit}.h5,t=${TI[i]}:${TF[i]}"
  if (( ${#ThinningPerFile[i]} )); then FitList="${FitList}t${ThinningPerFile[i]}"; fi
  FitList="${FitList},e=$NumExp$FitOptionsModel"
  if (( ${#FitOptionsPerFile[i]} )); then FitList="$FitList,${FitOptionsPerFile[i]}"; fi
  Title="$Title${MySep}ΔT=${DeltaT[i]}"
  LabelTI="$LabelTI${MySep}${PlotDataFrom}"
  LabelTF="$LabelTF${MySep}$(( ${DeltaT[i]} - PlotDataFrom ))"
  FitType=${FitType}_${TI[i]}_${TF[i]}
  MySep=" "
done

BuildModel=$OutSubDir/${FitWhat}_${OutLongName}

MultiFit="MultiFit -e 2 --Hotelling 0"
MultiFit="$MultiFit --overwrite"
[ -v FitOptions ] && MultiFit="$MultiFit $FitOptions"
[ -v UnCorr ] && MultiFit="$MultiFit --uncorr"
MultiFit="$MultiFit --debug-signals"
Cmd="$MultiFit --summary 2 -o $BuildModel $InputMesonSnk $InputMesonSrc $FitList"
BuildModel=$BuildModel.$FitType.$SuffixModel
#echo "A: $Cmd"
echo "$Cmd"  > $BuildModel.$MLUSeed.log
if  ! $Cmd &>> $BuildModel.$MLUSeed.log
then
  LastError=${PIPESTATUS[0]}
  if [ "$LastError" = 3 ]; then
    echo "Warning: Not all parameters resolved"
  else
    echo "Warning $LastError: $Cmd"
  fi
fi

# Get the fit characteristics: energy difference, matrix element, test stat, ...
Partial=EDiff,MEL${Gamma}0
[[ $FitWhat = R3 ]] && Partial=$Partial,R3${Gamma}Raw
ColumnValues=$(GetColumn --exact ChiSqPerDof,pValueH --partial $Partial ${BuildModel}.$MLUSeed.h5)
if [ "$?" != 0 ]; then
  echo "Error: $ColumnValues"
  unset ColumnValues
  unset RefText
else
  #echo "OK: $ColumnValues"
  ColumnValues=($ColumnValues)
  EDiff="${ColumnValues[@]:16:8}"
  MEL="${ColumnValues[@]:24:8}"
  [[ $FitWhat = R3 ]] && R3Raw="${ColumnValues[@]:32:8}"
  RefText="MEL${Gamma}0=${ColumnValues[25]} ${UnCorr+uncorrelated }χ²/dof=${ColumnValues[4]} (pH=${ColumnValues[12]})"
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

Cmd="ti='$LabelTI' tf='$LabelTF'"
[ -v yrangeFit ] && Cmd="$Cmd yrange='$yrangeFit'"
[ -v RefText ] && Cmd="$Cmd RefText='Fit $FitName: $RefText'"
[ -v RefValFit ] && Cmd="$Cmd RefVal='$RefValFit'"
Cmd="$Cmd title='$Title' field=$Field plottd.sh ${BuildModel}_td.$MLUSeed.txt"
#echo "B: $Cmd"
echo "$Cmd"  >> $BuildModel.$MLUSeed.log
eval  $Cmd  &>> $BuildModel.$MLUSeed.log

# Use the model to create the alternate

OutSubDirName=$OutSubDir/2Plot

unset MySep
unset FitList
for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  ThisFile=${PrefixPlot}${DeltaT[i]}${SuffixPlot}
  FitList="$FitList${MySep}${ThisFile}.h5,t=${TI[i]}:${TF[i]}${PlotOptions}"
  MySep=" "
done

mkdir -p $OutSubDirName
PlotModelBase=$OutSubDirName/${PlotWhat}_${OutLongName}

Cmd="$MultiFit -o $PlotModelBase ${BuildModel}.$MLUSeed.h5 $FitList"
PlotModelBase=$PlotModelBase.${FitType}
PlotModel=$PlotModelBase.$SuffixModel
#echo "C: $Cmd"
echo "$Cmd"  > $PlotModel.$MLUSeed.log
if !  $Cmd &>> $PlotModel.$MLUSeed.log
then
  LastError=${PIPESTATUS[0]}
  echo "Warning $LastError: $Cmd"
fi

# Plot it

Cmd="ti='$LabelTI' tf='$LabelTF'"
[ -v yrangePlot ] && Cmd="$Cmd yrange='$yrangePlot'"
[ -v RefText ] && Cmd="$Cmd RefText='Plot $PlotName: $RefText'"
[ -v RefValPlot ] && Cmd="$Cmd RefVal='$RefValPlot'"
Cmd="$Cmd title='$Title' field=$PlotField plottd.sh ${PlotModel}_td.$MLUSeed.txt"
#echo "D: $Cmd"
echo "$Cmd"      > $PlotModelBase.log
if ! eval $Cmd &>> $PlotModelBase.log
then
  LastError=${PIPESTATUS[0]}
  echo "Error $LastError: $Cmd"
fi
)

############################################################

# Perform simul fit

############################################################

# Parameters:
#   1: pSnk
# DeltaT: array of DeltaT to fit
# TI: array of fit start times
# TF: array of fit start times
# Thin: thinning or other options
# qSnk qSrc qSpec Gamma
# IncludeSpatial
# UnCorr: set to anything to perform uncorrelated fit
# FitOptions: extra options for MultiFit command-line
function DoSimulFit()
(
  local pSnk=$1
  local pSrc=0
  local MesonSnk; GetMesonFile MesonSnk $qSnk $qSpec
  local MesonSrc; GetMesonFile MesonSrc $qSrc $qSpec
  local i Title
  local FitBase="$Ensemble/MELFit/3sm_sp2/R3_${qSnk}_${qSrc}_${Gamma}_p2_${pSnk}.dt"
  if [[ $pSnk = 0 || $Gamma != gT ]]; then unset IncludeSpatial; fi
  local NumDeltaT=${#DeltaT[@]}
  mkdir -p "${FitBase%/*}"
  for(( i=0;i<NumDeltaT;++i)); do
    FitBase="${FitBase}_${DeltaT[i]}"
    Title="$Title${Title+ }'R3 ΔT=${DeltaT[i]}'"
  done
  FitBase="${FitBase}.Simul"
  local FitSuffix="${UnCorr+un}corr_${aKaonTIP[pSnk]}_${aKaonTFP[pSnk]}_${aKaonTIW[pSnk]}_${aKaonTFW[pSnk]}_${aDsTIP[pSrc]}_${aDsTFP[pSrc]}_${aDsTIW[pSrc]}_${aDsTFW[pSrc]}"
  local Files="$PlotData/corr/2ptp2/${MesonSnk}_p2_${pSnk}_g5P_g5P.fold.$MLUSeed.h5,t=${aKaonTIP[pSnk]}:${aKaonTFP[pSnk]},e=2"
    Files="${Files} $PlotData/corr/2ptp2/${MesonSnk}_p2_${pSnk}_g5P_g5W.fold.$MLUSeed.h5,t=${aKaonTIW[pSnk]}:${aKaonTFW[pSnk]},e=1"
    Files="${Files} $PlotData/corr/2ptp2/${MesonSrc}_p2_${pSrc}_g5P_g5P.fold.$MLUSeed.h5,t=${aDsTIP[pSrc]}:${aDsTFP[pSrc]},e=2"
    Files="${Files} $PlotData/corr/2ptp2/${MesonSrc}_p2_${pSrc}_g5P_g5W.fold.$MLUSeed.h5,t=${aDsTIW[pSrc]}:${aDsTFW[pSrc]},e=1"
  for(( i=0;i<NumDeltaT;++i)); do
    Files="${Files} $PlotData/ratioE1ZV1/3sm_sp2/R3_${qSnk}_${qSrc}_${Gamma}_dt_${DeltaT[i]}_p2_${pSnk}_g5P_g5P.fold.$MLUSeed.h5,t=${TI[i]}:${TF[i]}${Thin[i]},C2e=2,C2Model=cosh,raw"
    FitSuffix="${FitSuffix}_${TI[i]}_${TF[i]}"
  done
  if [ -v IncludeSpatial ]; then
    for(( i=0;i<NumDeltaT;++i)); do
      Files="${Files} $PlotData/ratioE1ZV1/3sm_sp2/R3_${qSnk}_${qSrc}_gXYZ_dt_${DeltaT[i]}_p2_${pSnk}_g5P_g5P.fold.$MLUSeed.h5,t=${TI[i]}:${TF[i]}${Thin[i]},C2e=2,C2Model=cosh,raw"
      FitSuffix="${FitSuffix}_${TI[i]}_${TF[i]}"
    done
  fi
  local OutBase="$FitBase.$FitSuffix.g5P_g5W.model"
  local LogFile="$OutBase.$MLUSeed.log"
  local Cmd="MultiFit -e $NumExp --Hotelling 0 --overwrite --debug-signals --summary 2"
  [ -v FitOptions ] && Cmd="$Cmd $FitOptions"
  [ -v UnCorr ] && Cmd="$Cmd --uncorr"
  Cmd="$Cmd -o $FitBase $Files"
  echo $Cmd &>  $LogFile
  if ! eval $Cmd &>> $LogFile; then echo "Returned ${PIPESTATUS[0]}" &>> $LogFile; fi

  # Get the fit characteristics: energy difference, matrix element, test stat, ...
  Exact=ChiSqPerDof,pValueH,${MesonSnk}_p2_${pSnk}-E0,${MesonSrc}_p2_${pSrc}-E0
  Partial=EDiff,R3${Gamma}Raw,MEL${Gamma}0
  if [ -v IncludeSpatial ]; then Partial=$Partial,R3gXYZRaw,MELgXYZ0; fi
  ColumnValues=$(GetColumn --exact $Exact --partial $Partial ${OutBase}.$MLUSeed.h5)
  if [ "$?" != 0 ]; then
    echo "Error: $ColumnValues"
    unset ColumnValues
    unset RefText
  else
    #echo "OK: $ColumnValues"
    ColumnValues=($ColumnValues)
    EDiff="${ColumnValues[@]:32:8}"
    R3Raw="${ColumnValues[@]:40:8}"
    MEL="${ColumnValues[@]:48:8}"
    Stats="χ²/dof=${ColumnValues[4]} (pH=${ColumnValues[12]})"
    #RefText="MEL${Gamma}0=${ColumnValues[25]} $Stats"
  fi

  local PlotFile="${OutBase}_td.$MLUSeed.txt"
  RefText="Fit R_3: K m_{eff}(${pSnk})=${ColumnValues[17]} $Stats" RefVal="${ColumnValues[@]:16:8}" title="'K p-p' 'K p-w'" yrange=${ayRange[$MesonSnk,$pSnk]} mmax=1 plottd.sh $PlotFile &>> $LogFile
  RefText="Fit R_3: D_s m_{eff}(${pSnk})=${ColumnValues[25]} $Stats" RefVal="${ColumnValues[@]:24:8}" title="'D_s p-p' 'D_s p-w'" yrange=${ayRange[$MesonSrc,$pSrc]} mmin=2 mmax=3 plottd.sh $PlotFile &>> $LogFile
  RefText="Fit R_3: MEL${Gamma}0(${pSnk})=${ColumnValues[49]} $Stats" RefVal="${ColumnValues[@]:40:8}" title="$Title" field=corr yrange=${ayrangeR3[$Gamma,$pSnk]} mmin=4 mmax=$((3+NumDeltaT)) plottd.sh $PlotFile &>> $LogFile
  if [ -v IncludeSpatial ]; then
    RefText="Fit R_3: MELgXYZ0(${pSnk})=${ColumnValues[65]} $Stats" RefVal="${ColumnValues[@]:56:8}" title="$Title" field=corr yrange=${ayrangeR3[gXYZ,$pSnk]} mmin=$((4+NumDeltaT)) mmax=$((3+2*NumDeltaT)) plottd.sh $PlotFile &>> $LogFile
  fi
)
