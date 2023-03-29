#!/bin/env bash

# Make an analysis job which will perform fits
# Optional environment variables:
# fixe:   Don't use an E fit file - instead set E=1
# fixzv:  Don't use a ZV fit file - instead set ZV=1
JobRoot=${JobRoot:-analyse}
InCorr=${InCorr:-corr}
InFit=${InFit:-fit}
InZVFit=${InZVFit:-ZVFit}
Out=${Out:-ratio}

# How to reach correlators during run
CorrDir=../$InCorr
FitDir=../$InFit
FitDirZV=../$InZVFit
# How to reach directories now (while this script runs)
CorrDirNow=$JobRoot/$InCorr
FitDirNow=$JobRoot/$InFit
FitDirZVNow=$JobRoot/$InZVFit
OutDirNow=$JobRoot/$Out

# Don't expect these to be overriden
JobBaseName=${OutDirNow//\//.}${fixe+E1}${fixzv+ZV1}
JobFileName=${JobBaseName}.sh

# Delete output file
[ -w $JobFileName    ] && rm $JobFileName

#set -x           #Log all steps - debug
set -e            #Fail on any errors
shopt -s nullglob #Return nothing if no match

function DoesFileExist()
{
    if ! [ -f "$1" ]
    then
	echo "$1 doesn't exist"
	unset bOK
    fi
}

# Walk each of the 3pt directories
for d in $CorrDirNow/3??_*
do
  Dir=${d##*/}
  Spec=${Dir:4}
  SpecGrp=$Spec
  unset PGroup
  if [ "${Spec: -2}" == "p2" ]
  then
    PGroup=${Spec: -2}
    Spec=${Spec:0:-2}
  fi
  # Work out whether we have the energy fits and the ZV Fits
  bOK=1
  if [ -v fixe ]
  then
    EFit="''" # ,${Spec}
  else
    EFit="Fit_${SpecGrp}.txt"
    DoesFileExist $FitDirNow/$EFit
    EFit=$FitDir/$EFit
  fi
  if [ -v fixzv ]
  then
    unset ZVFit
  else
    ZVFit="ZVFit_${SpecGrp}.txt"
    DoesFileExist $FitDirZVNow/$ZVFit
    ZVFit=$FitDirZV/$ZVFit
  fi
  if ! [ -v bOK ]
  then
    echo "Can't make ratios for $Dir"
  else
	ZVFit=R${ZVFit:+,$ZVFit}
	echo "Making ratios for $Dir"
	# Get a list (associative array) of the work to be done
	unset Meson
	unset DeltaT
	unset PW
	declare -A Meson
	declare -A PW
	for f in $d/*_gT_dt_*_p2_0_ps2_1_*.h5
	do
	    g=${f##*/}
	    Suffix=${g#*.}
	    g=${g%%.*}
	    parts=(${g//_/ })
	    if (( ${#parts[@]} == 12 )) && [ "${parts[1]}" != "${parts[2]}" ]
	    then
		Meson[${parts[1]}_${parts[2]}]=1
		PW[${parts[10]}_${parts[11]}]=1
	    fi
	done
	# Now make a job for each Quark
	Cmd="CRatio --i3 $CorrDir/$Dir/ --i2 $CorrDir/2pt${PGroup}/ -o $Dir/"
	Cmd="$Cmd --efit $EFit --type $ZVFit"
	unset f
	for m in ${!Meson[@]}
	do
	    for ss in ${!PW[@]}
	    do
		echo "$Cmd '*'_${m}_gT_dt_'*'_p2_0_ps2_'*'_${ss}.'*'.h5" >> $JobFileName
		echo "$Cmd '*'_${m}_gXYZ_dt_'*'_p2_0_ps2_'*'_${ss}.'*'.h5" >> $JobFileName
	    done
	done
  fi
done

# Job file needs to be executable
chmod u+x $JobFileName
