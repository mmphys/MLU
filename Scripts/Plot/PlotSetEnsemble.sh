#!/usr/bin/env bash

# Common routines used by plotting
set -e
#set -x

Ensemble="$1"
PlotDocDir=~/Work/Uni/PhD/SemiLep/Plot

# Get full path of parameter 2 in variable referenced by 1
function PCGetFullPath()
{
  local Dir
  local DirLen
  if [[ "$2" =~ .*/ ]]
  then
      Dir="${BASH_REMATCH[0]}"
      DirLen=${#Dir}
      if (( DirLen > 1 )) && [ ${Dir: -1} = "/" ]; then Dir="${Dir:0:$((DirLen-1))}"
      fi
  else
    Dir="."
  fi
  export $1="$( cd "$Dir" && pwd )"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  RetExit=exit # Being executed
else
  RetExit=return # Being sourced
fi

# Get the full path to this script
PCGetFullPath PlotCommonRoot "${BASH_SOURCE[0]}"
EnsembleBase=Ensemble

if ! [ -f "$PlotCommonRoot/$EnsembleBase$Ensemble.sh" ]; then
  Echo "$PlotCommonRoot/$EnsembleBase$Ensemble.sh missing"
  $RetExit 1
fi

if ! [ -f "$PlotDocDir/$EnsembleBase$Ensemble.tex" ]; then
  Echo "$PlotDocDir/$EnsembleBase$Ensemble.tex missing"
  $RetExit 1
fi

(
  cd $PlotCommonRoot
  [ -h "$EnsembleBase.sh" ] && rm "$EnsembleBase.sh"
  ln -s "$EnsembleBase$Ensemble.sh" "$EnsembleBase.sh"
)

(
  cd $PlotDocDir
  [ -h "$EnsembleBase.tex" ] && rm "$EnsembleBase.tex"
  ln -s "$EnsembleBase$Ensemble.tex" "$EnsembleBase.tex"
)
