#!/usr/bin/env bash

#set -x
set -e

MLUDir="$(cd "${0%/*}/.."; echo $PWD)"

function GetXcodeGUIBuild
{
  local SaveGlob="$(shopt -p globstar nullglob)"; shopt -s globstar nullglob
  local -n BuildDir=$1 # Returns XCode build directory for this program
  unset BuildDir
  local Package="${2:-$1}"
  local DerivedData=$HOME/Library/Developer/Xcode/DerivedData
  local BuildProd=Build/Products
  local -a Dirs=("$DerivedData"/"${Package}-"**/"$BuildProd")
  case ${#Dirs[@]} in
    0) echo "$Package no Xcode $BuildProd";;
    1) BuildDir="${Dirs[0]}"; echo "$Package in $BuildDir";;
    *) echo "Multiple $Package: ${Dirs[*]}";;
  esac
  eval "$SaveGlob"
}

# Not used – this returns command-line build location and I want GUI location
function GetBuildDir()
{
  local -n ActiveConfig=$1 # Returns the active config (usually Debug or Release)
  local -n BuildDir=$2 # Returns XCode build directory for this config
  local Dir="$3"
  local Project="$4"
  local XCProj="$5"
  local ProjPath="$Dir/$Project/$XCProj"
  local A=$(xcodebuild -project "$ProjPath" -showBuildSettings 2>/dev/null | grep TARGET_BUILD_DIR)
  BuildDir=${A#*= }
  if [ -z "$BuildDir" ]; then
    echo "Unable to find TARGET_BUILD_DIR for $ProjPath"
    #exit 1
    ActiveConfig=
    return 1
  else
    ActiveConfig="${BuildDir##*/}"
    #echo "Linking $GridPre$ActiveConfig to $Project $ActiveConfig in $BuildDir"
    echo "Linking $GridPre$ActiveConfig to $Project $ActiveConfig"
  fi
}

# Arg 1 is the config to install for (usually Debug or Release)
function CDInstallDir()
{
  local InstallDir="${GridPre:-.local}$1"
  mkdir -p "$InstallDir"
  cd "$InstallDir"
  mkdir -p bin lib include
}

function UnLink()
{
  local File
  for File in "$@"
  do
    if [ -e "$File" ] || [ -h "$File" ]; then rm "$File"; fi
  done
}

function LinkTo()
{
  local File Dir="$1"; shift
  if [ -n "$Dir" ]; then
    for File in "$@"
    do
      if [ -r "$Dir/$File" ]; then
        if [ -e "$File" ] || [ -h "$File" ]; then
          echo multiple copies of "$File"
        else
          ln -s "$Dir/$File"
        fi
      fi
    done
  fi
}

function DoConfig()
{
  local Config GridBP MLU SemiLep CGrid CMLU CSemiLep
  local GridLibs=(libGrid.a libHadrons.a)
  local MLULibs=(libMLU.dylib libMLU.0.dylib libMLU.a libMLU.la)
  local HadronsExe=(BatchDeflationBenchmark HadronsContractor \
            HadronsContractorBenchmark HadronsXmlRun HadronsXmlValidate \
            HadronsFermionEP64To32)
  local MLUExe=(bootstrap Continuum corr corrstat CRatio Divide FitSummary GetColumn ImportCorr \
            MakeEnsembleInfo mixedop MultiFit repair SynthData)
  local SemiLepExe=(Debug GaugeCmp hlxml peramb-compare TensorDemo xml3pt \
            xmlGFWall xmlPeramb xmlZ2)

  GetXcodeGUIBuild GridBP Grid
  GetXcodeGUIBuild MLU
  GetXcodeGUIBuild SemiLep
  for Config in "$@"
  do
    echo "$Config"
    #if GetBuildDir GridXConfig GridXBuild "$MikeSource" GridX Grid.xcodeproj; then
    CDInstallDir "$Config"
    [ -n "$GridBP" ] && CGrid="$GridBP/$Config"
    [ -n "$MLU" ] && CMLU="$MLU/$Config"
    [ -n "$SemiLep" ] && CSemiLep="$SemiLep/$Config"

    cd lib
    UnLink "${GridLibs[@]}" "${MLULibs[@]}"
    LinkTo "$CGrid" "${GridLibs[@]}"
    LinkTo "$CMLU" "${MLULibs[@]}"
    LinkTo "$CSemiLep" "${MLULibs[@]}"

    cd ../bin
    UnLink {grid,hadrons}-config "${HadronsExe[@]}" "${MLUExe[@]}" "${SemiLepExe[@]}"
    [ -n "$Grid" ] && LinkTo "$Grid/build$Config" grid-config
    [ -n "$Hadrons" ] && LinkTo "$Hadrons/build$Config" hadrons-config
    LinkTo "$CGrid" "${HadronsExe[@]}"
    LinkTo "$CMLU" "${MLUExe[@]}"
    LinkTo "$CSemiLep" "${MLUExe[@]}" "${SemiLepExe[@]}"

    cd ../include
    if [ -h MLU ]; then rm MLU
    elif [ -d MLU ]; then rm -rf MLU
    fi
    ln -s "$MLUDir/MLU"
    if [ -n "$Hadrons" ]; then
      if [ -h Hadrons ]; then rm Hadrons
      elif [ -d Hadrons ]; then rm -rf Hadrons
      fi
      ln -s "$Hadrons/Hadrons"
    fi
    # Grid
    if [ -n "$Grid" ]; then
      rm -rf Grid
      mkdir -p Grid
      cd Grid
      LinkTo "$Grid/build$Config/Grid" Config.h Version.h
      for f in "$Grid/Grid/"*; do if [[ -d $f ]] ; then rm -rf ${f##*/}; ln -s $f; fi; done
      for f in "$Grid/Grid/"*.{h,hpp}; do if [[ -f $f ]] ; then rm -f ${f##*/}; ln -s $f; fi; done
      cd ..
    fi
  done
}

if (($#))
  then DoConfig "$@"
else
  DoConfig Debug Release
fi
