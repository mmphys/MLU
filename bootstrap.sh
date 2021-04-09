#!/usr/bin/env bash

function Download() {
  local Location="$1"
  local Name="$2"
  local Zip="$3"
  if [ "$Zip" == "" ]; then Zip="${Name}.tar.gz"; fi

  echo "... Installing $Name"
  if ! wget $Location/$Zip
  then
    # CERN math libraries appear to be offline. Fallback to my copy
    wget https://www2.ph.ed.ac.uk/~s1786208/$Zip
  fi
  local val=$?
  if [ $? -ne 0 ]
  then
    echo "Couldn't download $Name from $Location/$Zip"
  else
    tar -xf $Zip
    echo "$Name installed"
  fi
  return $val
}

# Fail on any error
set -e

echo "Meson Lattice Utilities (MLU)"
if [[ -d .Package ]]
then
  echo "Skipping pre-requisites (delete PreReq to reinstall)"
else
  echo "Installing Pre-requisites"
  mkdir .Package
  cd .Package
  if Download https://mirror.ibcp.fr/pub/gnu/gsl "GSL v2.6" gsl-2.6.tar.gz \
    && Download http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2 "Minuit2 v5.34.14" Minuit2-5.34.14.tar.gz
  then
    echo "Prerequisites installed"
  else
    echo "Prerequisites not installed"
  fi

  cd ..
fi

echo "Generating configure for Meson Lattice Utilities (MLU)"
autoreconf -fvi
