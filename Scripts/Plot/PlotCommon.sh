#!/usr/bin/env bash

# Common routines used by plotting
set -e
#set -x

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  echo "${BASH_SOURCE[0]} should be sourced"
  set +e; exit 1
fi

# Copy the specified file to the same directory with a fixed name
# Parameters:
#   1: Full path to file to copy
#   2: Standardised name for file, default A_Latest
function PCSaveLatest()
{
  set -x
  local InFile=$1
  local StandardName=${2-A_Latest}

  local Ext=${InFile##*/}; Ext=${Ext##*.}
  local OutFile=${InFile%/*}; OutFile=$OutFile/$StandardName.$Ext
  echo cp $InFile $OutFile
  set +x
}

# Get full path
# Parameters:
#   1: Filename to get full path of
#   2: Variable to return value in (default FullPath)
function PCGetFullPath()
{
  declare -n RetVal=${2:-FullPath}
  local Dir
  local DirLen
  if [[ "$1" =~ .*/ ]]
  then
      Dir="${BASH_REMATCH[0]}"
      DirLen=${#Dir}
      if (( DirLen > 1 )) && [ ${Dir: -1} = "/" ]; then Dir="${Dir:0:$((DirLen-1))}"
      fi
  else
    Dir="."
  fi
  RetVal="$( cd "$Dir" && pwd )"
}

# Get quark name without mass suffix
# Parameters:
#   1: Quark to get name of
#   2: Variable to return value in (default qName)
function PCQNameNoMass()
{
  declare -n RetVal=${2:-qName}
  local q="$1"
  case "${q:0:1}" in
    l) q=light;;
    s) q=strange;;
    h) q=heavy;;
  esac
  RetVal="$q"
}

# Get quark name with mass suffix
# Parameters:
#   1: Quark to get name of
#   2: Variable to return value in (default qName)
function PCQName()
{
  declare -n RetVal=${2:-qName}
  PCQNameNoMass "$1" "${!RetVal}"
  if ! [ -z "${1:1}" ]; then
    RetVal="$RetVal (m=0.${1:1})"
  fi
}

# Get human readable version of sink and source
# Parameters:
#   1: Operator at sink (e.g. g5P)
#   2: Operator at source (e.g. g5P)
#   3: Variable to return value in (default opHuman)
function PCPointWall()
{
  declare -n RetVal=${3:-opHuman}
  local opSnk="$1"
  local opSrc="$2"
  case ${opSnk}_${opSrc} in
    g5P_g5P) RetVal="point-point";;
    g5P_g5W) RetVal="point-wall";;
    g5W_g5P) RetVal="wall-point";;
    g5W_g5W) RetVal="wall-wall";;
          *) RetVal="${opSnk}-${opSrc}";;
  esac
}

# Get column values from .h5 file
# Parameters:
#   1: Filename
#   2: Prefix for RefText
#   3: Comma separated list of exact field names
#   4: Comma separated list of partial field names
# Optional:
#   UnCorr: Set to anything to indicate uncorrelated fit
# Returns:
#   ColumnValues: Array of column values in following order:
#                 1) ChiSqPerDof 2) pValueH 3) Exact parameters 4) Partial parameters
#   RefText: Reference string containing first column value, chi^2/dof and pValueH
#   On error, these are unset
function GetColumnValues()
{
  local Filename="$1"
        RefText="$2"
  local Exact="$3"
  local Partial="$4"
  if ColumnValues=$(GetColumn --exact "ChiSqPerDof,pValueH${Exact:+,$Exact}" \
                    ${Partial:+--partial "$Partial"} "$Filename")
  then
    #echo "OK: $ColumnValues"
    ColumnValues=($ColumnValues)
    RefText="$RefText${ColumnValues[@]:17:1},"
    [ -v UnCorr ] && RefText="$RefText uncorrelated"
    RefText="$RefText χ²/dof=${ColumnValues[@]:4:1}, pH=${ColumnValues[@]:12:1}"
  else
    LastError=${PIPESTATUS[0]}
    #echo "Error $LastError: $ColumnValues" # GetColumn shows its own error messages
    unset ColumnValues
    unset RefText
    #return $LastError
  fi
}

# Get script name without path prefix or extension
function PCMakeScriptOutDir()
{
  local OutDir="${BASH_SOURCE[1]}"
  OutDir="${OutDir##*/}"
  OutDir="${OutDir%.*}"
  OutDir="${OutDir#Plot}"
  mkdir -p "$Ensemble/$OutDir"
  cd "$Ensemble/$OutDir"
}

# File name already split apart
function HumanReadable
{
  DeltaTHuman="ΔT=${DeltaT}"
  case $Gamma in
      gT) GammaHuman="V_0";;
    gXYZ) GammaHuman="V_i";;
       *) GammaHuman="$Gamma";;
  esac
  PCPointWall ${opSnk} ${opSrc}
}

# 1: Where to save the Meson (file) name
# 2: q
# 3: spectator
function GetMesonFile()
{
  declare -n MesonFile=$1
  local q="$2"
  local Spec="$3"
  case "${q:0:1}_${Spec:0:1}" in
    s_h | l_h | l_s) MesonFile="${Spec}_${q}";;
                  *) MesonFile="${q}_${Spec}";;
  esac
}

# 1: Where to save the Meson (human-readable) name
# 2: q
# 3: spectator
function GetMeson()
{
  declare -n NameHuman=$1
  local q="$2"
  local Spec="$3"
  local qNum
  [ -v OptionNoMass ] && unset qNum || qNum="${q:1}"
  case "${q:0:1}_${Spec:0:1}" in
    h_s | s_h) NameHuman="D_s";;
    h_l | l_h) NameHuman="D";;
    l_l)       [ -v OptionASCII ] && NameHuman="pi" || NameHuman="π";;
    l_s | s_l) NameHuman="K";;
            *) NameHuman="${q}-${Spec}";;
  esac
  if ! [ -v OptionASCII ] && ! [ -z "$qNum" ]; then NameHuman="$NameHuman ($qNum)"; fi
}

# 1: 3pt filename to split
# 2: Spectator
function Split3ptFile()
{
  local f="${1##*/}"
  local Spec="$2"
  Ext=${f#*.}
  Prefix=${f%%.*}
  PrefixParts=(${Prefix//_/ })
  local NumParts=${#PrefixParts[@]}
  if (( NumParts != 12 && NumParts != 10 )) || [[ "${PrefixParts[4]}" != dt ||
    "${PrefixParts[6]}" != p2 || ( $NumParts = 12 && "${PrefixParts[8]}" != ps2 ) ]]
  then
    #echo "${#PrefixParts[@]}: ${PrefixParts[4]}, ${PrefixParts[6]}, ${PrefixParts[8]}."
    echo "Bad 3pt filename $f"
    return 1
  fi
  Ratio=${PrefixParts[0]} #quark or anti
  qSnk=${PrefixParts[1]}
  qSrc=${PrefixParts[2]}
  Gamma=${PrefixParts[3]}
  DeltaT=${PrefixParts[5]}
  p2=${PrefixParts[7]}
  (( NumParts == 12 )) && ps2=${PrefixParts[9]} || unset ps2
  opSnk=${PrefixParts[-2]}
  opSrc=${PrefixParts[-1]}
  HumanReadable
  GetMeson MSnkHuman $qSnk $Spec
  GetMeson MSrcHuman $qSrc $Spec
  pMax=$p2
  if (( pMax < ps2 )); then pMax=$ps2; fi
}

# Return the longer of two strings
function PCGetLonger()
{
  local A="$1"
  local B="$2"
  declare -n RetVal=${3:-Longer}
  if (( ${#B} > ${#A} )); then
    RetVal="$B"
  else
    RetVal="$A"
  fi
}

############################################################

# Get the full path to this script
PCGetFullPath "${BASH_SOURCE[0]}" PlotCommonRoot

# Load parameters
PlotEnsemble="$PlotCommonRoot/Ensemble${Ensemble}.sh"
if ! [ -e "$PlotEnsemble" ]; then
  echo "Missing: $PlotEnsemble"
  set +e; return 1
fi
unset PlotData
. "$PlotEnsemble"
if ! [ -v PlotData ]; then
  echo "$PlotEnsemble did not define PlotData"
  set +e; exit 1
elif ! [ -d "$PlotData" ]; then
  echo "Ensemble $Ensemble: $PlotData missing"
  set +e; exit 1
fi

# Computed parameters
LHalf=$((L/2))
THalf=$((T/2))
export MLUSeed=${MLUSeed:-1835672416} # Seed to process data
DataSeed=${DataSeed:-$MLUSeed} # Seed in data filenames
