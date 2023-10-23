#!/usr/bin/env bash

# Functions for different types of two-point fits (source me)

export size="${size:-5.8in,3in}"
export sizeSimulP="${sizeSimulP:-5.8in,4in}"

############################################################

# Execute command

# Input
#  Cmd:     Command to execute
#  LogFile: Where to redirect output

############################################################

function DoCmd()
{
  if [[ "$LogFile" != "${LogFile%/*}" && -n "${LogFile%/*}" ]]; then mkdir -p "${LogFile%/*}"; fi
  if [ "$1" = EraseLog ] && [ -f $LogFile ]; then rm $LogFile; fi
  #echo $Cmd
  echo $Cmd &>> $LogFile
  local RetVal=0
  if ! ( eval $Cmd ) &>> $LogFile; then
    RetVal=${PIPESTATUS[0]}
    # Ignore MultiFit 'Not all params resolved" warning
    if [[ ${Cmd:0:8} = MultiFit && $RetVal = 3 ]]; then
      RetVal=0
      echo "Warning: not all parameters resolved: $Cmd"
    else
      echo "Error $RetVal" &>> $LogFile
      echo "Error $RetVal executing $Cmd"
    fi
  fi
  return $RetVal
}

############################################################

# Perform a simple, two-point fit of point-point and point-wall

# Optional:
#   ExtraName: Series name added to end of fit base name
#   Stat:      Minimum acceptable pvalueH

############################################################

function FitTwoPoint()
{
########################################
# Input parameters
local qHeavy=${qHeavy:-h$Heavy}
local Meson=${Meson:-${qHeavy}_s}
local p=${p:-0}
local NumExp=${NumExp:-3}
local NumExp2=${NumExp2:-${NumExp}}
local MaxExp=$(( NumExp2 > NumExp ? NumExp2 : NumExp ))
local TI=${TI:-4}
local TF=${TF:-18}
local TI2=${TI2:-$TI}
local TF2=${TF2:-$TF}
local MaxTF=$(( TF2 > TF ? TF2 : TF ))
# Stat - which statistic
# FitOptions - extra options to MultiFit
local Thin1=${Thin1+,t$Thin1}
local Thin2=${Thin2+,t$Thin2}

#LabelTF=${LabelTF:-$((THalf-2))}
local LabelTF=${LabelTF:-$(( MaxTF + 2 ))}
# LabelTI also used (optional)

local Corr=${Corr:-corr}

########################################
# Derived

local OutSubDir=$Ensemble/MELFit/2ptp2/$Meson
local CorrPrefix=$Corr/2ptp2/${Meson}_p2_${p}_
local Suffix=fold.$DataSeed

mkdir -p "$OutSubDir"

########################################
# Perform the fit

local Cmd="MultiFit -e $NumExp --summary 2"
[ -v Stat ] && Cmd+=" --Hotelling $Stat"
[ -v ExtraName ] && Cmd+=" --extra '$ExtraName'"
Cmd+=" --overwrite"
Cmd+=" --debug-signals"
Cmd+=" --strict 1"
[ -v FitOptions ] && Cmd+=" $FitOptions"
Cmd+=" -i $PlotData/${CorrPrefix} -o $OutSubDir/"
Cmd+=" g5P_g5P.$Suffix.h5,t=${TI}:${TF}${Thin1},e=$NumExp"
if (( TF2 >= TI2 )); then
  Cmd+=" g5P_g5W.$Suffix.h5,t=${TI2}:${TF2}${Thin2},e=$NumExp2"
fi

#echo "$Cmd"
local ModelBase="$(eval $Cmd --showname)"
local LogFile=$ModelBase.$MLUSeed.log
local FitFile=$ModelBase.$MLUSeed.h5
local TDFile=${ModelBase}_td.$MLUSeed.txt
#echo "--showname=$ModelBase"

if [ -v PlotOnly ]; then
  LogFile=${LogFile%.log}.plot.log
  echo "Skipping fit: $Cmd" > $LogFile
else
DoCmd EraseLog || return 0
fi

########################################
# Plot the fit

# Get E0 (reference value for plot)
local MesonHuman
OptionNoMass= GetMeson MesonHuman ${Meson//_/ }
Latex= GetColumnValues $FitFile '$'"$MesonHuman"'(n^2='"$p"'), \\, E_0=$' '' E0
Cmd="title='point-point point-wall' tf='$LabelTF $LabelTF' yrange='${ayRange[$Meson,$p]}'"
[ -v LabelTI ] && Cmd+=" ti='$LabelTI $LabelTI'"
[ -v RefText ] && Cmd+=" RefText='$RefText' RefVal='${ColumnValues[@]:16:8}'"
local ExtraFiles="$PlotData/${CorrPrefix}gT5P_g5P.$Suffix.txt"
ExtraFiles+=" $PlotData/${CorrPrefix}gT5P_g5W.$Suffix.txt"
Cmd+=" extra='$ExtraFiles'"
Cmd+=' Latex= ylabel='"'"'$a E_\\textrm{eff}$'"'"
Cmd+=" plottd.sh $TDFile"
DoCmd
}

############################################################

# Perform a fit of multiple correlators
# Parameters
#   FileList: Array of files to fit
#   FileArgs: Array of fit arguments for each file
#      See GetColumnValues for more detail on params 1, 2 & 3

# Optional:
#   *:         See plottd.sh for all variables it supports
#   FileTitle: Array of titles for each file
#   ExtraName: Series name added to end of fit base name
#   Stat:      Minimum acceptable pvalueH
#   NumExp:    Default number of exponentials
#   FitOptions: Extra options for MultiFit
#   RefText:   Prefix text for reference value (default E_0=)
#   RefField:  Which field to return as reference (default E0)

# Returns:
#   ColumnValues: Array of reference value fields

############################################################

function FitMulti()
{
########################################
# Input parameters
local RefText="${RefText-E_0=}"
local InPrefix="${InPrefix-$PlotData/corr/}"
local OutPrefix="$Ensemble/$OutPrefix"
local Suffix="${Suffix-.fold.$DataSeed}"
local RefField="${RefField:-E0}"

if ! (( ${#FileList[@]} )); then echo "FitMulti() FileList empty"; return; fi

########################################
# Derived

#mkdir -p "$OutPrefix"

########################################
# Generic fit command
local Cmd="MultiFit"
Cmd+=" --summary 2"
Cmd+=" --overwrite"
Cmd+=" --debug-signals"
Cmd+=" --strict 1"
[ -v Stat ] && Cmd+=" --Hotelling $Stat"
[ -v ExtraName ] && Cmd+=" --extra '$ExtraName'"
[ -v NumExp ] && Cmd+=" -e $NumExp"
[ -v FitOptions ] && Cmd+=" $FitOptions"
[ -n "$InPrefix" ] && Cmd+=" -i '$InPrefix'"
[ -n "$OutPrefix" ] && Cmd+=" -o '$OutPrefix'"

# Build list of files to fit
local title
for((i=0; i<${#FileList[@]}; ++i))
do
  Cmd+=" '${FileList[i]}$Suffix.h5'"
  (( ${#FileArgs[i]} )) && Cmd+=",${FileArgs[i]}"
  # Plot title for this file
  (( ${#title} )) && title+=' '
  (( ${#FileTitle[i]} )) && title+="'${FileTitle[i]}'" || title+="'${FileList[i]}'"
done

# Get the name of the resulting fit file
#echo "$Cmd"
local ModelBase="$(eval $Cmd --showname)"
local LogFile=$ModelBase.$MLUSeed.log
local FitFile=$ModelBase.$MLUSeed.h5
local TDFile=${ModelBase}_td.$MLUSeed.txt
#echo "--showname=$ModelBase"

# Perform the fit
if [ -v PlotOnly ]; then
  LogFile=${LogFile%.log}.plot.log
  echo "Skipping fit: $Cmd" > $LogFile
else
DoCmd EraseLog || return 0
fi

########################################
# Plot the fit

# Get reference value for plot
GetColumnValues $FitFile "$RefText" '' "$RefField"
Cmd=". plottd.sh $TDFile"
DoCmd
}

############################################################

# Scan a two-point fit range

############################################################

function TwoPointScan()
{
  local OptP="$1" # t=...[,e=...]
  local OptW="$2"

  local qHeavy=${qHeavy:-h$Heavy}
  local Meson=${Meson:-${qHeavy}_s}
  local p=${p:-0}

  # SubDir - which subdirectory to place output in
  # SummaryOptions - extra options for the Summary

  local Corr=${Corr:-corr}

  ########################################
  # Error check

  local DoErrorMsg=0
  if [ -z "$OptP" ]
  then
    echo "Point options unspecified. Expect as a minimum t=ti:tf"
    DoErrorMsg=1
  fi
  if [ -z "$OptW" ]
  then
    echo "Wall options unspecified. Expect as a minimum t=ti:tf"
    DoErrorMsg=1
  fi
  if (( DoErrorMsg )); then return 1; fi

  ########################################
  # Derived

  if [[ ! $p =~ _ ]]; then p="p2_$p"; fi # Assume momentum is p^2 (but allow it to be set explicitly)

  local InDir=$PlotData/$Corr/2ptp2
  local OutDir=$Ensemble/TwoPointScan/$Meson
  if [ -v SubDir ]; then OutDir="$OutDir/$SubDir"; fi

  local Suffix=fold.$DataSeed.h5
  local FileBase=${Meson}_${p}_g5P_g5
  local LogFile=$OutDir/${FileBase}W.log

  ########################################
  # Perform the fit

  echo "$Ensemble: Scan two point ${Meson}_${p} $SubDir, Point $OptP, Wall $OptW"

  mkdir -p $OutDir # So I can redirect output here

  local Cmd="MultiFit -o \"$OutDir/\" -i \"$InDir/$FileBase\""
  Cmd="$Cmd --debug-signals --mindp 2 --iter 100000 --strict 3"
  [ -v FitOptions ] && Cmd="$Cmd $FitOptions"
  Cmd="$Cmd P.$Suffix,$OptP"
  Cmd="$Cmd W.$Suffix,$OptW"
  DoCmd || return 0

  Cmd="FitSummary --all --strict 0 -i \"$OutDir/\" -o \"$OutDir/\"" # --errdig 2 --tablen
  [ -v SummaryOptions ] && Cmd="$Cmd $SummaryOptions"
  Cmd="$Cmd '${Meson}_${p}.*.h5'"
  DoCmd

  Cmd="rm '$OutDir/${Meson}_${p}.corr_'*"
  DoCmd
}

############################################################

# Scan standard ranges around start and stop times for point and wall

############################################################

function StdTwoPointScan()
{
  local TIP=$1
  local TFP=$2
  local TIW=$3
  local TFW=$4
  local NumExp=${NumExp:-2}
  local NumExp2=${NumExp2:-1}
  SubDir=${SubDir}a TwoPointScan t=$((TIP-2)):$((TFP-1)):6:3,e=${NumExp} t=${TIW}:${TFW},e=${NumExp2}
  SubDir=${SubDir}b TwoPointScan t=${TIP}:${TFP},e=${NumExp} t=$((TIW-2)):$((TFW-1)):6:3,e=${NumExp2}
}


############################################################

# Simultaneous fits of point-point data at all momenta

############################################################

function SimulP()
{
  local Meson=${1:-s_l}
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local i EKeys title
  sTimes="${aTimes[*]}"
  sTimes=${sTimes// /_}
  sTimes=${sTimes//:/_}
  Cmd="$MultiFit -e 2 -N $L -o $OutDir/$Meson"
  [ -n "$ExtraName" ] && Cmd+=" --extra '$ExtraName'"
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do
    Cmd+=" ${aFitFiles[i]}.h5,t=${aTimes[i]}${aThin[i]:+t${aThin[i]}}"
  done
  
  #echo "$Cmd"
  local ModelBase="$(eval $Cmd --showname)"
  local LogFile=$ModelBase.$MLUSeed.log
  local FitFile=$ModelBase.$MLUSeed.h5
  local TDFile=${ModelBase}_td.$MLUSeed.txt
  #echo "--showname=$ModelBase"

  if [ -v PlotOnly ]; then
    LogFile=${LogFile%.log}.plot.log
    echo "Skipping fit: $Cmd" > $LogFile
  else
  DoCmd EraseLog || return 0
  fi

  # Get the energy difference
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do
    EKeys+=${EKeys:+,}${Meson}_p2_${i}-E0
  done
  local MesonHuman
  OptionNoMass= GetMeson MesonHuman ${Meson//_/ }
  Latex= GetColumnValues $FitFile '$'"$MesonHuman"'(n^2=0), \\, E_0=$' $EKeys

  # Now plot it
  for((i=0; i < ${#aFitFiles[@]}; ++i)); do title+="${title:+ }"'$n^2='"$i"'$'; done
  Cmd="title='$title'"
  [ -v RefText ] && Cmd+=" RefText='$RefText'"
  [ -v sizeSimulP ] && Cmd+=" size='$sizeSimulP'"
  Cmd+=' Latex= ylabel='"'"'$a E_\\textrm{eff}$'"'"
  Cmd+=" plottd.sh $TDFile"
  DoCmd

  local ModelBaseName=${ModelBase##*/}
  ModelBaseName=${ModelBaseName#*.}
  for((i=0; i < ${#aFitFiles[@]}; ++i))
  do
    Cmd="title='n^2=$i'"
    [ -v RefText ] && Cmd+=" RefText='E_0(n^2=${i})=${ColumnValues[@]:$((17+i*8)):1}, χ²/dof=${ColumnValues[@]:4:1} (pH=${ColumnValues[@]:12:1})'"
    Cmd+=" yrange='${ayRange[s_l,$i]}'"
    Cmd+=" save='$OutDir/${Meson}_p2_$i.${ModelBaseName}_log.$MLUSeed'"
    Cmd+=" mmin=$i mmax=$i plottd.sh $TDFile"
    DoCmd
  done
}

############################################################

# For each non-zero momentum
# fit pp & pw data using E(n^2=0) and lattice dispersion relation

# Input:
#  aTimes:  array of point-point start:stop fit times for each momenta (0'th ignored)
#  aTimesW: array of point-wall  start:stop fit times for each momenta (0'th ignored)
#  PriorFitTimes: fit to use for zero momentum

############################################################

function FitEachMomPW()
{
  local OutGroup=$1
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local ZeroMomFit="$OutDir/s_l_p2_0.corr_$PriorFitTimes.g5P_g5W.model.1835672416.h5"
  for((i=1; i < 5; ++i))
  do
    sTimes="${aTimes[i]} ${aTimesW[i]}"
    sTimes=${sTimes// /_}
    sTimes=${sTimes//:/_}
    BaseFile="$OutDir/s_l_p2_$i.${OutGroup}_$PriorFitTimes"
    Cmd="$MultiFit -e 2 -N 24 -o $BaseFile $ZeroMomFit"
    ThisFile="${aFitFiles[i]}.h5"
    Cmd="$Cmd $ThisFile,t=${aTimes[i]}"
    Cmd="$Cmd ${ThisFile//g5P_g5P/g5P_g5W},t=${aTimesW[i]},e=1"
    BaseFile="$BaseFile.corr_$sTimes.g5P_g5W.model"
    LogFile="$BaseFile.$MLUSeed.log"
    if [ -v PlotOnly ]; then
      LogFile=${LogFile%.log}.plot.log
      echo "Skipping fit: $Cmd" > $LogFile
    else
    DoCmd || return 0
    fi

    GetColumnValues $BaseFile.$MLUSeed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    Cmd="title='point-point point-wall' "
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$MLUSeed.txt"
    DoCmd
  done
}

############################################################

# For each non-zero momentum
# fit pp & pw data using E(n^2=0) and lattice dispersion relation

# Input:
#  aTimes:  array of point-point start:stop fit times for each momenta (0'th ignored)
#  PriorFitTimes: fit to use for zero momentum

############################################################

function FitEachMomP()
{
  local OutGroup=$1
  local Cmd
  local sTimes
  local BaseFile
  local LogFile
  local ZeroMomFit="$OutDir/s_l_p2_0.corr_$PriorFitTimes.g5P_g5W.model.1835672416.h5"
  for((i=1; i < 5; ++i))
  do
    sTimes="${aTimes[i]}"
    sTimes=${sTimes//:/_}
    BaseFile="$OutDir/s_l_p2_$i.${OutGroup}_$PriorFitTimes"
    Cmd="$MultiFit -e 2 -N 24 -o $BaseFile $ZeroMomFit"
    Cmd="$Cmd ${aFitFiles[i]}.h5,t=${aTimes[i]}"
    BaseFile="$BaseFile.corr_$sTimes.g5P_g5W.model"
    LogFile="$BaseFile.$MLUSeed.log"
    if [ -v PlotOnly ]; then
      LogFile=${LogFile%.log}.plot.log
      echo "Skipping fit: $Cmd" > $LogFile
    else
    DoCmd || return 0
    fi

    GetColumnValues $BaseFile.$MLUSeed.h5 "E_0(n^2=$i)=" s_l_p2_${i}-E0

    Cmd="title='point-point n^2=$i'"
    [ -v RefText ] && Cmd="$Cmd RefText='$RefText'"
    Cmd="$Cmd yrange='${ayRange[s_l,$i]}'"
    Cmd="$Cmd plottd.sh ${BaseFile}_td.$MLUSeed.txt"
    DoCmd
  done
}
