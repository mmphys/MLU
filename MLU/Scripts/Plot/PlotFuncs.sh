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
  local -a dT=(${BaseParts[@]:3})
  local i Meson ti='' tf='' title='' DeltaTOna='\\flatfrac{\\Delta T}{a}'
  OptionNoMass= GetMeson Meson $Q $Spec
  [ "$Meson" == π ] && Meson='\\pi'

  # Get the fit characteristics: energy difference, matrix element, test stat, ...
  local Latex=
  GetColumnValues "$FitFile" "\$Z_V($Meson)=\$" '' ZV

  # tMin / TMax ranges
  for (( i=0; i<${#dT[@]}; ++i )); do
    if ((i)); then
      ti+=' '
      tf+=' '
    fi
    ti+='2'
    tf+="$((dT[i]-2))"
  done

  # Title
  #if ((${#dT[@]}>1)); then
    for (( i=0; i<${#dT[@]}; ++i )); do
      ((i)) && title+=' '
      title+="'\$$DeltaTOna=${dT[i]}\$'"
    done
  #elif ! [ -v MikeThesis ]; then
    #title+="'"
    #title+='$Z_V \\left[ '"$Meson, "' \\,'
    #[ "$Meson" == K ] && title+=' \\textrm{spec=}'"$Spec, "' \\,'
    #title+="$DeltaTOna=${dT[0]}, "' \\,'
    #title+='\\textrm{'"${PW//_/\\\\_}"'} '
    #title+=' \\right] $'
    #title+="'"
  #fi
  if [ -n "$title" ]; then
    export title
    echo "title=$title"
  fi

  # Plot it
  
  Cmd='field=corr'
  [ -v yrange ] && Cmd+=" yrange='$yrange'"
  if [ -v RefText ]; then
    Cmd+=" RefText='$RefText' RefVal='${ColumnValues[@]:16:8}'"
  fi
  Cmd+=" ti='$ti' tf='$tf'"
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
