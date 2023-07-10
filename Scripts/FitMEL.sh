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
# FitOptionsRatio (optional) extra fit options for ratios, e.g. C2eSrc=3
# FileSeries
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

DeltaT=(${DeltaT:-24 28 32})
TI=(${TI:-13 14 14})
TF=(${TF:-16 20 24})

case "$FitWhat" in
  R3 | quark) : ;;
          "") FitWhat=R3;;
           *) echo "FitWhat=$FitWhat unrecognised"; exit 1;;
esac
NumExp=${NumExp:-3} # How many exponentials in the 3pt fit
FileSeries=${Alt+alt_}${FileSeries:-E$NumExp}

PW=${PW:-g5P_g5P}
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

local PlotOptions FitOptionsModel

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

local FileOpModel; PCGetLonger $FileOpSnk $FileOpSrc FileOpModel
SuffixModel=$FileOpModel.model
MesonDir=$Ensemble/$MELFit/2ptp2
InputMesonSnk=$MesonDir/$MesonSnk/${MesonSnk}${FileMomSnk}.$FitSnk.${FileOpSnk}.model.$MLUSeed.h5
InputMesonSrc=$MesonDir/$MesonSrc/${MesonSrc}${FileMomSrc}.$FitSrc.${FileOpSrc}.model.$MLUSeed.h5

local OutBaseName=${qSnk}_${qSrc}_${Gamma}_p2_$((pSrc > pSnk ? pSrc : pSnk))
local ExtraName=dt_$DeltaTAll.$FitSnk.$FitSrc${FileSeries+.$FileSeries}
local OutSubDir=$Ensemble/$MELFit/$SpecDir

############################################################

# Main loop

############################################################

# Fit the data

mkdir -p $OutSubDir

local MySep
for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  FitList+="$MySep$PrefixFit${DeltaT[i]}$SuffixFit.h5,e=$NumExp$FitOptionsModel,t=${TI[i]}:${TF[i]}"
  [ -n "${ThinningPerFile[i]}" ] && FitList+="t${ThinningPerFile[i]}"
  [ -n "${FitOptionsPerFile[i]}" ] && FitList+=",${FitOptionsPerFile[i]}"
  [ -n "$FitOptionsRatio" ] && FitList+=",$FitOptionsRatio"
  Title+="${MySep}ΔT=${DeltaT[i]}"
  LabelTI+="$MySep${PlotDataFrom}"
  LabelTF+="$MySep$(( ${DeltaT[i]} - PlotDataFrom ))"
  MySep=" "
done

MultiFit="MultiFit -e 2 --Hotelling 0 --overwrite --debug-signals --summary 2 --extra '$ExtraName'"
[ -v UnCorr ] && MultiFit+=" --uncorr"
[ -v FitOptions ] && MultiFit+=" $FitOptions"
Cmd="$MultiFit -o $OutSubDir/${FitWhat}_$OutBaseName $InputMesonSnk $InputMesonSrc $FitList"

#echo "A: $Cmd"
local ModelBase="$(eval $Cmd --showname)"
local LogFile=$ModelBase.$MLUSeed.log
local FitFile=$ModelBase.$MLUSeed.h5
local TDFile=${ModelBase}_td.$MLUSeed.txt

     echo "$Cmd"  > $LogFile
if ! eval  $Cmd &>> $LogFile
then
  LastError=${PIPESTATUS[0]}
  if [ "$LastError" = 3 ]; then
    echo "Warning: Not all parameters resolved"
  else
    echo "Warning $LastError: $Cmd"
  fi
fi

# Get the fit characteristics: energy difference, matrix element, test stat, ...
local RefValFit RefValPlot
Partial=EDiff,MEL${Gamma}0
[[ $FitWhat = R3 ]] && Partial=$Partial,R3${Gamma}Raw
ColumnValues=$(GetColumn --exact ChiSqPerDof,pValueH --partial $Partial $FitFile)
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
  case "$FitWhat" in
    quark)
          RefValFit="$EDiff";;
    R3)
          RefValFit="$R3Raw"
          RefValPlot="$EDiff";;
  esac
fi

# Plot it
Cmd="ti='$LabelTI' tf='$LabelTF'"
[ -v yrangeFit ] && Cmd+=" yrange='$yrangeFit'"
[ -v RefText ] && Cmd+=" RefText='Fit $FitName: $RefText'"
[ -v RefValFit ] && Cmd+=" RefVal='$RefValFit'"
Cmd+=" title='$Title' field=$Field plottd.sh $TDFile"
#echo "B: $Cmd"
echo "$Cmd"  >> $LogFile
eval  $Cmd  &>> $LogFile

# Use the model to create the alternate

unset MySep
unset FitList
for (( i = 0; i < ${#DeltaT[@]}; ++i ))
do
  FitList+="$MySep$PrefixPlot${DeltaT[i]}$SuffixPlot.h5,e=$NumExp$PlotOptions,t=${TI[i]}:${TF[i]}"
  [ -n "${ThinningPerFile[i]}" ] && FitList+="t${ThinningPerFile[i]}"
  MySep=" "
done

OutSubDirPlot=$OutSubDir/2Plot
mkdir -p $OutSubDirPlot

Cmd="$MultiFit -o $OutSubDirPlot/${PlotWhat}_$OutBaseName $FitFile $FitList"
PlotModel=$PlotModelBase.$SuffixModel

#echo "C: $Cmd"
local PlotBase="$(eval $Cmd --showname)"
      LogFile=$PlotBase.$MLUSeed.log
local PlotTDFile=${PlotBase}_td.$MLUSeed.txt

     echo "$Cmd"  > $LogFile
if ! eval  $Cmd &>> $LogFile
then
  LastError=${PIPESTATUS[0]}
  echo "Warning $LastError: $Cmd"
fi

# Plot it

Cmd="ti='$LabelTI' tf='$LabelTF'"
[ -v yrangePlot ] && Cmd+=" yrange='$yrangePlot'"
[ -v RefText ] && Cmd+=" RefText='Plot $PlotName: $RefText'"
[ -v RefValPlot ] && Cmd+=" RefVal='$RefValPlot'"
Cmd+=" title='$Title' field=$PlotField plottd.sh $PlotTDFile"
#echo "D: $Cmd"
echo "$Cmd"     >> $LogFile
if ! eval $Cmd &>> $LogFile
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
# SourcePriorFit: set to prior Source fit
# SinkPriorDisp: set to anything to use the dispersion relation for sink
function DoSimulFit()
(
  local pSnk=$1
  local pSrc=0
  local MesonSnk; GetMesonFile MesonSnk $qSnk $qSpec
  local MesonSrc; GetMesonFile MesonSrc $qSrc $qSpec
  local i Title
  local FitBaseDir="$Ensemble/MELFit/3sm_sp2"
  local FitBase="$FitBaseDir/R3_${qSnk}_${qSrc}_${Gamma}_p2_${pSnk}"
  if [[ $pSnk = 0 || $Gamma != gT ]]; then unset IncludeSpatial; fi
  local NumDeltaT=${#DeltaT[@]}
  if [ -v SinkPriorDisp ] && ! [ -v SourcePriorFit ]; then
    echo "Can't use dispersion for sink without a prior fit"
    return 1
  fi
  mkdir -p "$FitBaseDir"
  local ExtraName=dt
  for(( i=0;i<NumDeltaT;++i)); do
    ExtraName+="_${DeltaT[i]}"
    Title+="${Title+ }'R3 ΔT=${DeltaT[i]}'"
  done
  ExtraName+=".Simul"
  local Files
  if ! [ -v SinkPriorDisp ]; then
    Files+=" $PlotData/corr/2ptp2/${MesonSnk}_p2_${pSnk}_g5P_g5P.fold.$MLUSeed.h5,t=${aKaonTIP[pSnk]}:${aKaonTFP[pSnk]},e=2"
    Files+=" $PlotData/corr/2ptp2/${MesonSnk}_p2_${pSnk}_g5P_g5W.fold.$MLUSeed.h5,t=${aKaonTIW[pSnk]}:${aKaonTFW[pSnk]},e=1"
  fi
  if [ -v SourcePriorFit ]; then
    Files+=" $FitBaseDir/$SourcePriorFit,${MesonSrc}_p2_0-E,${MesonSrc}_p2_0-g5P"
    if [ -v SinkPriorDisp ]; then
      # TODO: This next only works for Z2 overlap coeffficients which are momentum independent
      Files+=",${MesonSnk}_p2_0-E,${MesonSnk}_p2_0-g5P=${MesonSnk}_p2_${pSnk}-g5P"
    fi
  else
    Files+=" $PlotData/corr/2ptp2/${MesonSrc}_p2_${pSrc}_g5P_g5P.fold.$MLUSeed.h5,t=${aDsTIP[pSrc]}:${aDsTFP[pSrc]},e=2"
    Files+=" $PlotData/corr/2ptp2/${MesonSrc}_p2_${pSrc}_g5P_g5W.fold.$MLUSeed.h5,t=${aDsTIW[pSrc]}:${aDsTFW[pSrc]},e=1"
  fi
  for(( i=0;i<NumDeltaT;++i)); do
    Files+=" $PlotData/ratioE1ZV1/3sm_sp2/R3_${qSnk}_${qSrc}_${Gamma}_dt_${DeltaT[i]}_p2_${pSnk}_g5P_g5P.fold.$MLUSeed.h5,t=${TI[i]}:${TF[i]}${Thin[i]},C2e=2,C2Model=cosh,raw${FitOptionsRatio+,$FitOptionsRatio}"
  done
  if [ -v IncludeSpatial ]; then
    for(( i=0;i<NumDeltaT;++i)); do
      Files+=" $PlotData/ratioE1ZV1/3sm_sp2/R3_${qSnk}_${qSrc}_gXYZ_dt_${DeltaT[i]}_p2_${pSnk}_g5P_g5P.fold.$MLUSeed.h5,t=${TI[i]}:${TF[i]}${Thin[i]},C2e=2,C2Model=cosh,raw${FitOptionsRatio+,$FitOptionsRatio}"
    done
  fi
  local Cmd="MultiFit -e $NumExp --Hotelling 0 --overwrite --debug-signals --summary 2"
  [ -v SinkPriorDisp ] && Cmd+=" -N $L"
  [ -v FitOptions ] && Cmd+=" $FitOptions"
  [ -v UnCorr ] && Cmd+=" --uncorr"
  Cmd+=" -o $FitBase --extra '$ExtraName' $Files"

  #echo $Cmd
  local ModelBase="$(eval $Cmd --showname)"
  local LogFile=$ModelBase.$MLUSeed.log
  local FitFile=$ModelBase.$MLUSeed.h5
  local TDFile=${ModelBase}_td.$MLUSeed.txt

  echo $Cmd &>  $LogFile
  if ! eval $Cmd &>> $LogFile; then echo "Returned ${PIPESTATUS[0]}" &>> $LogFile; fi

  # Get the fit characteristics: energy difference, matrix element, test stat, ...
  Exact=ChiSqPerDof,pValueH,${MesonSnk}_p2_${pSnk}-E0,${MesonSrc}_p2_${pSrc}-E0
  Partial=EDiff,R3${Gamma}Raw,MEL${Gamma}0
  if [ -v IncludeSpatial ]; then Partial=$Partial,R3gXYZRaw,MELgXYZ0; fi
  ColumnValues=$(GetColumn --exact $Exact --partial $Partial $FitFile)
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

  ModelDTStart=0
  if ! [ -v SinkPriorDisp ]; then
    RefText="Fit R_3: K m_{eff}(${pSnk})=${ColumnValues[17]} $Stats" RefVal="${ColumnValues[@]:16:8}" title="'K p-p' 'K p-w'" yrange=${ayRange[$MesonSnk,$pSnk]} mmax=1 plottd.sh $TDFile &>> $LogFile
    ModelDTStart=$((ModelDTStart+2))
  fi
  if ! [ -v SourcePriorFit ]; then
    RefText="Fit R_3: D_s m_{eff}(${pSnk})=${ColumnValues[25]} $Stats" RefVal="${ColumnValues[@]:24:8}" title="'D_s p-p' 'D_s p-w'" yrange=${ayRange[$MesonSrc,$pSrc]} mmin=$ModelDTStart mmax=$((ModelDTStart+1)) plottd.sh $TDFile &>> $LogFile
    ModelDTStart=$((ModelDTStart+2))
  fi
  RefText="Fit R_3: MEL${Gamma}0(${pSnk})=${ColumnValues[49]} $Stats" RefVal="${ColumnValues[@]:40:8}" title="$Title" field=corr yrange=${ayrangeR3[$Gamma,$pSnk]} mmin=$ModelDTStart mmax=$((ModelDTStart+NumDeltaT-1)) plottd.sh $TDFile &>> $LogFile
  if [ -v IncludeSpatial ]; then
    RefText="Fit R_3: MELgXYZ0(${pSnk})=${ColumnValues[65]} $Stats" RefVal="${ColumnValues[@]:56:8}" title="$Title" field=corr yrange=${ayrangeR3[gXYZ,$pSnk]} mmin=$((ModelDTStart+NumDeltaT)) mmax=$((ModelDTStart+2*NumDeltaT-1)) plottd.sh $TDFile &>> $LogFile
  fi
)
