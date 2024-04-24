#!/usr/bin/env bash

# User input
# Grid:        Path to Grid
# Hadrons:     Path to Hadrons
# GridPre:     Install prefix
# Headers:     Set to anything to install links to headers only

if [ -z "$GridPre" ]; then echo '$GridPre not set'; exit 1; fi

# Make wildcard failures return empty strings (not the unmodified wildcard string)
shopt -s nullglob

function LinkFiles()
{
  local Src Dst
  local Dir="$1"; shift
  mkdir -p "$Dir"
  cd "$Dir"
  for Src in "$@"
  do
    Dst="${Src##*/}"
    rm -f "$Dst"
    ln -s "$MyEnv/$Src"
  done
}

if [ -z "$Grid" ]; then
  echo '$Grid not set'
else
  echo 'Grid source located in '$Grid
  for MyEnv in $Grid/build${GridSelect}*
  do
    ShortEnv=${MyEnv#${Grid}/build}
    Dest=$GridPre$ShortEnv
    echo "  ${ShortEnv}: $Dest --> $MyEnv"
    if ! [ -v Headers ]; then
      LinkFiles "$Dest/bin" grid-config
      LinkFiles "$Dest/lib" "Grid/libGrid.a"
    fi
    # Make Grid subdirectory with all required headers
    LinkFiles "$Dest/include/Grid" Grid/Config.h Grid/Version.h
    for f in $Grid/Grid/*; do if [[ -d $f ]] ; then rm -rf ${f##*/}; ln -s $f; fi; done
    for f in $Grid/Grid/*.{h,hpp}; do if [[ -f $f ]] ; then rm -f ${f##*/}; ln -s $f; fi; done
  done
fi

if [ -z "$Hadrons" ]; then
  echo '$Hadrons not set'
else
  echo 'Hadrons source located in '$Hadrons
  for MyEnv in $Hadrons/build${GridSelect}*
  do
    ShortEnv=${MyEnv#${Hadrons}/build}
    Dest=$GridPre$ShortEnv
    echo "  ${ShortEnv}: $Dest --> $MyEnv"
    # bin
    if ! [ -v Headers ]; then
      LinkFiles "$Dest/bin" hadrons-config utilities/Hadrons{Contractor,XmlRun,XmlValidate}
      LinkFiles "$Dest/lib" Hadrons/libHadrons.a
    fi
    LinkFiles "$Dest/include"
    rm -rf Hadrons
    ln -s $Hadrons/Hadrons
  done
fi
