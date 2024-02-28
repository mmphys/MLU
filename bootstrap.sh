#!/usr/bin/env bash

function Download() {
  local Name="$1"
  local Location="${2%/*}"
  local Zip="${2##*/}"
  local Sha256Value="$3"
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
      Location="http://lqcd.me/PhD/tar"
      echo "  Downloading $Zip from $Location"
      if ! wget $Location/$Zip; then return 1; fi
    fi
  fi

  # Check checksum
  if command -v sha256sum; then
   echo "$Sha256Value $Zip" | sha256sum --check || exit 1
  else
     echo "WARNING: could not verify checksum, please install sha256sum" >&2
  fi

  # Now unzip it
  echo "  $Name unzipping ..."
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
echo "Downloading Pre-requisites"
mkdir -p .Package
cd .Package
if  ! Download "GSL v2.7" \
      https://mirror.ibcp.fr/pub/gnu/gsl/gsl-2.7.tar.gz \
      efbbf3785da0e53038be7907500628b466152dbc3c173a87de1b5eba2e23602b \
 || ! Download "Minuit2 v5.34.14" \
      http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz \
      2ca9a283bbc315064c0a322bc4cb74c7e8fd51f9494f7856e5159d0a0aa8c356
# NB: Minuit2 standalone on cern seems to be password protected. But download defaults to my copy
# FYI Minuit2 source code here https://github.com/root-project/root/tree/master/math/minuit2
then
  echo "Prerequisites not installed" >&2
  exit 1
fi
cd ..
echo "Prerequisites installed"

echo "Generating configure for Meson Lattice Utilities (MLU)"
autoreconf -fvi
