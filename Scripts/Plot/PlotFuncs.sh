#!/usr/bin/env bash

############################################################

# Plot functions sourced from PlotCommon.sh

############################################################

############################################################

# Plot ZV fit

# 1: ZV fit file (.h5) to plot
# Spec: Spectator (optional)

############################################################

function PlotZVFit()
{
  local FitFile="$1"
  local Spec="$2"

  # Get spectator from filename if not caller provided
  local PGroup
  if [ -z "$Spec" ]; then
    GetSpectator "$FitFile"
    if [ -z "$Spec" ]; then
      echo "Unknown spectator $FitFile"
      return 1
    fi
  fi

  local TDFile="${FitFile//.model./.model_td.}"
        TDFile="${TDFile%.h5}.txt"
  local Filename="${FitFile##*/}"
  local FileParts=(${Filename//./ })
  local PW=${FileParts[2]}
  local Basename="${FileParts[0]}"
        Basename="${Basename%%.*}"
  local BaseParts=(${Basename//_/ })
  local Q=${BaseParts[1]}
  local dT=${BaseParts[3]}
  local Meson
  OptionNoMass= GetMeson Meson $Q $Spec
  [ "$Meson" == π ] && Meson='\\pi'

  # Get the fit characteristics: energy difference, matrix element, test stat, ...
  local Latex=
  GetColumnValues "$FitFile" "\$Z_V($Meson)=\$" '' ZV

  # Plot it
  
  Cmd="yrange='$yrange' field=corr"
  if ! [ -v MikeThesis ]; then
    Cmd+=' title="'"'"'\$Z_V \\\\left[ $Meson, '
    [ "$Meson" == K ] && Cmd+='\\\\textrm{ spec=}$Spec, '
    Cmd+='\\\\, \\\\Delta T=$dT, \\\\, \\\\textrm{$PW} \\\\right]\$'"'"'"'
  fi
  if [ -v RefText ]; then
    Cmd+=" RefText='$RefText' RefVal='${ColumnValues[@]:16:8}'"
  fi
  Cmd+=" ti=2 tf=$((dT-2))"
  Cmd+=" size='${size:-5in,1.8in}'"
  Cmd+=" Latex="
  if [ -v MikeThesis ]; then
    Cmd+=" ylabel='\$Z_V\$'"
  else
    Cmd+=" ylabel='\$Z_V = \\\\flatfrac{\\\\tilde{C}^{(2)}\\\\left(\\\\Delta T\\\\right)}{C^{(3)}\\\\left(t\\\\right)}\$'"
  fi
  Cmd+=" plottd.sh $TDFile"

  echo "$Cmd"
  eval  $Cmd
}