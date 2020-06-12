#!/usr/bin/env bash

if [[ "$Grid"    == "" ]] ; then echo '$Grid not set'   ; exit 1; fi
if [[ "$Hadrons" == "" ]] ; then echo '$Hadrons not set'; exit 1; fi
if [[ "$GridPre" == "" ]] ; then echo '$GridPre not set'; exit 1; fi

# Make wildcard failures return empty strings (not the unmodified wildcard string)
shopt -s nullglob

# Grid
echo 'Grid source located in '$Grid
for MyEnv in $Grid/build${GridSelect:-*}
do
  ShortEnv=${MyEnv#${Grid}/build}
  Dest=$GridPre/$ShortEnv
  echo "${ShortEnv}: $Dest --> $MyEnv"
  # bin
  mkdir -p $Dest/bin
  cd       $Dest/bin
  rm -f grid-config
  ln -s $MyEnv/grid-config
  # lib
  mkdir -p $Dest/lib
  cd       $Dest/lib
  rm -f libGrid.a
  ln -s $MyEnv/Grid/libGrid.a
  #ln -s $MyEnv/Hadrons/libHadrons.a
  # include
  for sub in Grid #Hadrons
  do
    mkdir -p $Dest/include/$sub
    cd       $Dest/include/$sub
    for f in $Grid/$sub/*; do if [[ -d $f ]] ; then rm -rf ${f##*/}; ln -s $f; fi; done
    for f in $Grid/$sub/*.{h,hpp}; do if [[ -f $f ]] ; then rm -f ${f##*/}; ln -s $f; fi; done
    if [[ $sub == "Grid" ]]
    then
      rm -f Config.h
      ln -s $MyEnv/$sub/Config.h
      rm -f Version.h
      ln -s $MyEnv/$sub/Version.h
    fi
  done
done

# Hadrons
echo 'Hadrons source located in '$Hadrons
for MyEnv in $Hadrons/build{GridSelect:-*}
do
  ShortEnv=${MyEnv#${Hadrons}/build}
  Dest=$GridPre/$ShortEnv
  echo "${ShortEnv}: $Dest --> $MyEnv"
  # bin
  mkdir -p $Dest/bin
  cd       $Dest/bin
  rm -f hadrons-config
  ln -s $MyEnv/hadrons-config
  for f in $MyEnv/utilities/Hadrons{Contractor,XmlRun,XmlValidate}
  do
    rm -f ${f##*/}
    ln -s $f
  done
  # lib
  mkdir -p $Dest/lib
  cd       $Dest/lib
  rm -f libHadrons.a
  ln -s $MyEnv/Hadrons/libHadrons.a
  # include
  mkdir -p $Dest/include
  cd       $Dest/include
  rm -rf Hadrons
  ln -s $Hadrons/Hadrons
done
