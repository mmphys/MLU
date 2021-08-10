#!/usr/bin/env bash

function Download() {
  local Name="$1"
  local Location="$2"
  local Zip="$3"
  #if [ "$Zip" == "" ]; then Zip="${Name}.tar.gz"; fi

  # Download the file from requested location, falling back to my copy if unavailable
  echo "Installing pre-requisite $Name"
  if [ -r "$Zip" ]; then
    echo "  $Zip already installed"
  else
    echo "  Downloading $Zip from $Location"
    if ! wget $Location/$Zip
    then
      # If unavailable from primary location, fallback to my copy
      Location="https://www2.ph.ed.ac.uk/~s1786208"
      echo "  Downloading $Zip from $Location"
      if ! wget $Location/$Zip; then return 1; fi
    fi
  fi

  # Now unzip it
  tar -xf $Zip
  local val=$?
  if [ $? -ne 0 ]
  then
    echo "  Couldn't unzip $Zip"
  else
    echo "  $Name unzipped"
  fi
  return $val
}

# Fail on any error
set -e

echo "Meson Lattice Utilities (MLU)"
if [[ -d .Package ]]
then
  echo "Warning: .Package (pre-requisite) directory already exists"
fi
  echo "Downloading Pre-requisites"
  mkdir -p .Package
  cd .Package
  if Download "GSL v2.7" https://mirror.ibcp.fr/pub/gnu/gsl gsl-2.7.tar.gz \
    && Download "Minuit2 v5.34.14" http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2 Minuit2-5.34.14.tar.gz
  then
    echo "Prerequisites installed"
  else
    echo "Prerequisites not installed"
  fi

  cd ..

echo "Generating configure for Meson Lattice Utilities (MLU)"
autoreconf -fvi
