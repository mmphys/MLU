#!/usr/bin/env bash

# Common routines used by plotting
set -e
#set -x

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  echo "${BASH_SOURCE[0]} should be sourced"
  set +e; exit 1
fi

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

# Get humand readable version of sink and source
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
            *) NameHuman="${q}${Spec}";;
  esac
  if ! [ -v OptionASCII ] && ! [ -z "$qNum" ]; then NameHuman="$NameHuman (m=0.$qNum)"; fi
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

# Get the full path to this script
PCGetFullPath "${BASH_SOURCE[0]}" PlotCommonRoot

# Load parameters
PlotEnsemble="$PlotCommonRoot/Ensemble.sh"
if ! [ -h "$PlotEnsemble" ]; then
  echo "Link missing: $PlotEnsemble"
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
