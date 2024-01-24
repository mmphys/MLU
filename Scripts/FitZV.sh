#!/usr/bin/env bash

############################################################

# Functions for performing various versions of matrix element fits

############################################################

. PlotCommon.sh

#set -x
set -e

############################################################

# Initialise variables from user

############################################################

if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit 2
fi
DefaultSeries=Const # Original ZV fits from constant
#DefaultSeries=Jan24 # Jan 2024 ZV Fits from model including excited-states
OutDir=$Ensemble/${OutDir:-Renorm}
ZVDir=${ZVDir:-ZV}
ZVPlotDir=${ZVPlotDir:-${ZVDir}Plot}
ZVEFit=E_For_${ZVDir}.txt
ZVLink=ZV.txt # Make this link to the default series

if [ -v DoAll ]; then
  DoEFit=
  DoZV=
  DoPlot=
  DoFit=
fi

############################################################

# Derived variables

############################################################

eval Actions=(h$Heavy s l)
eval Spectators=(l s h${Heavy}) # The heavy is a bit sh#t as a spectator

# y range for ZV for each action,spectator combination
declare -A aYRange

MultiFit="MultiFit --overwrite --debug-signals --strict 3 --summary 2"

############################################################

# Initialisation

############################################################

mkdir -p $OutDir
cd $OutDir

############################################################

# Perform ZV Energy fit

# Optional
# UnCorr:   uncorrelated fit
# Stat:     Hotelling p-value
# FitOptions: extra fit options

############################################################

function FitZVEnergy()
{
local Q1=$1
local Q2=$2
local yrange=$3
local ti=${ti:-6}
local tf=${tf:-28}
local ti2=${ti2:-$ti}
local tf2=${tf2:-$tf}
local e=${e:-2}
local e2=${e2:-$e}
local p2=0

# Derived
local Meson; OptionNoMass= GetMeson Meson ${Q1} ${Q2}
local InDir=$PlotData/corr/2ptp2/${Q1}_${Q2}_p2_${p2}_g5P_g5
local InSuffix=.fold.$DataSeed.h5
local OutDir=Fit/ZVEnergy

mkdir -p $OutDir

# Make Fit command
local Cmd="$MultiFit${Stat:+ --Hotelling $Stat}"
[ -v ExtraName ] && Cmd+=" --extra '$ExtraName'"
[ -v FitOptions ] && Cmd+=" $FitOptions"
[ -v UnCorr ] && Cmd+=" --uncorr"
Cmd+=" -o $OutDir/ -i $InDir P$InSuffix,t=${ti}:${tf},e=$e W$InSuffix,t=${ti2}:${tf2},e=${e2}"

#echo "A: $Cmd"
local ModelBase="$(eval $Cmd --showname)"
local LogFile=$ModelBase.$MLUSeed.log
local FitFile=$ModelBase.$MLUSeed.h5
local TDFile=${ModelBase}_td.$MLUSeed.txt
#echo "--showname=$ModelBase"

if [ -v PlotOnly ]; then
  LogFile=${LogFile%.log}.plot.log
  echo "Skipping fit: $Cmd" > $LogFile
else
# Execute Fit command
echo "$Cmd"  > $LogFile
if  ! $Cmd &>> $LogFile
then
  LastError=${PIPESTATUS[0]}
  if [ "$LastError" = 3 ]; then
    echo "Warning: Not all parameters resolved"
  else
    echo "Warning $LastError: $Cmd"
  fi
fi

# Add this to my fit list
[ -v ZVEFit ] && echo "${Q1}_${Q2}_p2_${p2} $FitFile" >> $ZVEFit
fi

# Get the fit characteristics: energy difference, matrix element, test stat, ...
Latex= GetColumnValues $FitFile "\\\$$Meson\\\$"' \$m_{\\textrm{eff}}=\$' '' E

# Plot it

Cmd="title=\"'\\\$$Meson\\\$ p-p' '\\\$$Meson\\\$ p-w'\""
Cmd+=" yrange='$yrange'"
if [ -v RefText ]; then
  Cmd+=" RefText='$RefText' RefVal='${ColumnValues[@]:16:8}'"
fi
Cmd+=" ti='$PlotTI' tf='$PlotTF'"
Cmd+=" size='${size:-5in,2in}'"
# ylabel start
Cmd+=" ylabel='"'\$a '
((p2)) && Cmd+='E' || Cmd+='m'
Cmd+='_{\\textrm{eff}}'
if ! [ -v MikeThesis ]; then
  Cmd+=' = \\ln\\left[ \\flatfrac{C(t-\\frac{1}{2})}{C(t+\\frac{1}{2})} \\right]'
fi
# ylabel end
Cmd+=' \$'"'"
Cmd+=" Latex="
Cmd+=" plottd.sh $TDFile"
#echo "B: $Cmd"
echo "$Cmd"  >> $LogFile
eval  $Cmd  &>> $LogFile
}

############################################################

# Perform ZV Energy fit

# Optional
# UnCorr:   uncorrelated fit

############################################################

function MakeZV()
{
local CRatio="CRatio --i2 $PlotData/corr/2ptp2/ --efit E_For_ZV.txt"
local Spectator Action InDir Cmd FileSpec
local LogFile=$ZVDir/MakeZV.log
mkdir -p $ZVDir
[ -e $LogFile ] && rm $LogFile
for Spectator in ${Spectators[@]}; do
  InDir=$PlotData/corr/3pt_${Spectator}p2
  Cmd="$CRatio --i3 $InDir/ -o $ZVDir/${Spectator}p2/"
  for Action in ${Actions[@]}; do
    FileSpec=*_${Action}_${Action}_gT_dt_*_p2_0_ps2_0_*.$DataSeed.h5
    if ls $InDir/$FileSpec &> /dev/null
    then
      echo "$Cmd $FileSpec" >> $LogFile
      eval  $Cmd $FileSpec &>> $LogFile
    fi
  done
done
}

############################################################

# Main loop

# Preferred ZV_strange uses light spectator, otherwise strange spectator

############################################################

function PlotZV()
{
local Edge=${Edge:-2}
local cond="column(1)<$Edge || column(1)>word(FileDT,File)-$Edge"
local x='column(1)-0.5*word(FileDT,File)+(word(FileDT,File)-24)*.025'
local fields=corr
local key='bottom center maxrows 2'
local yrange SpecDir dT legend FileName FileNames
export cond x fields key legend

mkdir -p $ZVPlotDir
for Action in ${Actions[@]}; do
  for Spectator in ${Spectators[@]}; do
    SpecDir=${Spectator}p2
    Combo=$Action,$Spectator
    [ -n "${aYRange[$Combo]}" ] && export yrange=${aYRange[$Combo]} || unset yrange
    for PW in g5{P,W}_g5{P,W}; do
      legend=
      FileNames=
      for dT in ${EnsembleDeltaT[@]}; do
        FileName=$ZVDir/$SpecDir/ZV_${Action}_dt_${dT}_${PW}.fold.$MLUSeed.txt
        if ls $FileName &> /dev/null; then
          legend="$legend${legend+ }ΔT=$dT"
          FileNames="$FileNames${FileNames+ }$FileName"
        fi
      done
      if [ -n "$FileNames" ]
      then
        save=$ZVPlotDir/ZV_${Action}_spec_${Spectator}_${PW} plot.sh $FileNames
      fi
    done
  done
done
}

############################################################

# Perform ZV fit
# Parameters - DeltaT list

# Mandatory:
#  Q:       Quark action we're getting Z_V for
#  ti:      Array of initial fit times
#  tf:      Array of final   fit times
# Optional
#  Spec:      Spectator we want to use (default 's')
#  Thin:      Array of thinning specs
#  Const:     Set to anything to perform a fit to constant
#  UnCorr:    uncorrelated fit
#  Stat:      Hotelling p-value
#  FitOptions: extra fit options
#  yrange:    gnuplot yrange
#  Snk:       Sink   (default: g5P)
#  Src:       Source (default: g5P)

############################################################

function FitZVNew()
{
local Spec=${Spec:-s}
local -a dT=($*) ti=($ti) tf=($tf)
local e=${e:-2}
local Snk=${Snk-g5P}
local Src=${Src-g5P}
if (( ${#dT[@]} < 1 || ${#dT[@]} != ${#ti[@]} || ${#dT[@]} != ${#tf[@]} )); then
  echo "${#dT[@]} deltaT but ${#ti[@]} ti and ${#tf[@]} tf"
  return 1
fi

# Derived
local i EFit ReadKey ReadVal ReadJunk ExtraName
local Meson; OptionNoMass= GetMeson Meson $Q $Spec
local PW=${Snk}; [ "$Snk" != "$Src" ] && PW=${PW}_${Src}
local SpecDir=${Spec}p2
local InDir=$ZVDir/$SpecDir
local OutDir=Fit/$ZVDir$Series/$SpecDir
local InPrefix=$InDir/ZV_${Q}_dt_
local InSuffix=_${Snk}_${Src}.fold.$MLUSeed.h5
local ZVName="$Q-ZV"

(( ${#dT[@]} > 1 )) && ExtraName=$(IFS=_; echo "_${dT[*]:1}")

# Unless fitting to constant, we need the energy fit
if ! [ -v Const ]; then
  while read ReadKey ReadVal ReadJunk; do
    if [ "$ReadKey" == "${Q}_${Spec}_p2_0" ] || [ "$ReadKey" == "${Spec}_${Q}_p2_0" ]
    then
      EFit="$ReadVal"
      break
    fi
  done < "$ZVEFit"
  if [ -z "$EFit" ]; then
    echo "Energy fit for ${Q}_$Spec not found in $ZVEFit"
    return 1
  fi
fi

mkdir -p $OutDir

# Make Fit command
local Cmd="$MultiFit"
[ -v Stat ] && Cmd+=" --Hotelling '$Stat'"
[ -n "$ExtraName" ] && Cmd+=" --extra '$ExtraName'"
[ -v FitOptions ] && Cmd+=" $FitOptions"
[ -v UnCorr ] && Cmd+=" --uncorr"
[ -v Const ] || Cmd+=" -e $e"
Cmd+=" -o $OutDir/"
[ -n "$EFit" ] && Cmd+=" $EFit"
for (( i=0; i<${#dT[@]}; ++i ))
do
  Cmd+=" $InPrefix${dT[i]}${InSuffix},t=${ti[i]}:${tf[i]}"
  [ -n "${Thin[i]}" ] && Cmd+="t${Thin[i]}"
  [ -v Const ] && Cmd+=",model=const,const=$ZVName" || Cmd+=',C2e=2'
done

#echo "A: $Cmd"
local ModelBase="$(eval $Cmd --showname)"
local LogFile=$ModelBase.$MLUSeed.log
local FitFile=$ModelBase.$MLUSeed.h5
local TDFile=${ModelBase}_td.$MLUSeed.txt
#echo "--showname=$ModelBase"

if [ -v PlotOnly ]; then
  LogFile=${LogFile%.log}.plot.log
  echo "Skipping fit: $Cmd" > $LogFile
else
# Execute Fit command
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
fi

# Add this to my fit list
[ -v ZVFit ] && echo "$Q"$'\t'"$FitFile" >> $ZVFit

#echo "B: $Cmd"
PlotZVFit "$FitFile" "$Spec" &>> $LogFile
}

############################################################

# Perform ZV fit to constant

# Optional
# UnCorr:   uncorrelated fit
# Stat:     Hotelling p-value
# FitOptions: extra fit options

############################################################

function FitZV()
{
local Q=$1
local dT=${2:-${EnsembleDeltaT[0]}}
local Spec=${3:-s}
# DEPRECATED - None of 4, 5 or 6 are ever used and ti/tf always set by caller
#local yrange=$4
#local Snk=${5:-g5P}
#local Src=${6:-g5P}
#local ti=${ti:-$((dT/2-1))}
#local tf=${tf:-$((dT/2+1))}
Const= FitZVNew $dT
}
