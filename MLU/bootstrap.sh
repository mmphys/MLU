#!/usr/bin/env bash

# Fail on any error
set -e

cd "${0%/*}"
echo "Bootstrapping Meson Lattice Utilities (MLU) in $PWD"

function Download() {
  local Name="$1"
  local Location="${2%/*}"
  local Zip="${2##*/}"
  local Sha256Value="$3"
  local Unzip="${4+N}"
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
   echo "$Sha256Value $Zip" | sha256sum --check || return 1
  else
     echo "WARNING: could not verify checksum, please install sha256sum" >&2
  fi

  # Now unzip it
  if [ -z "$Unzip" ]
  then
    echo "  $Name unzipping ..."
    if ! tar -xf $Zip
    then
      echo "  Couldn't unzip $Zip"
      return 2
    fi
    echo "  $Name unzipped"
  fi
}

echo "Downloading Pre-requisites"
AllOK=1
GSLPrefix=gsl-2.7
Minuit2Prefix=Minuit2-5.34.14
if  ! Download "GSL v2.7" \
      https://mirror.ibcp.fr/pub/gnu/gsl/${GSLPrefix}.tar.gz \
      efbbf3785da0e53038be7907500628b466152dbc3c173a87de1b5eba2e23602b
then
  AllOK=0
else
  cd "$GSLPrefix"
  echo "Configuring $GSLPrefix"
  autoreconf -fvi
  cd ..
  #if [ ! -d gsl ]
  #then
    #mkdir gsl
    #cd gsl
    #echo "Making links to gsl headers"
    #find ../gsl-2.7 -type f -name '*'.h -exec ln -s {} \;
    #cd ..
  #fi
  #[ ! -h gsl ] && ln -s "$GSLPrefix"/gsl
fi

# NB: Minuit2 standalone on cern seems to be password protected. But download defaults to my copy
# FYI Minuit2 source code here https://github.com/root-project/root/tree/master/math/minuit2
if  ! Download "Minuit2 v5.34.14" \
      http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/${Minuit2Prefix}.tar.gz \
      2ca9a283bbc315064c0a322bc4cb74c7e8fd51f9494f7856e5159d0a0aa8c356 \
      No
then
  AllOK=0
fi

if(( !AllOK ))
then
  echo "Prerequisites not installed" >&2
  exit 1
fi
echo "Prerequisites installed"

echo "Generating configure for Meson Lattice Utilities (MLU)"
autoreconf -fvi
