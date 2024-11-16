#!/usr/bin/env bash

# Make R1/R2/R3 and Z_V ratios with energy (&Z_V) frozen
RootDir=analyse
Corr=corr

if ! [ -d $RootDir/$Corr ]; then
  echo "Correlator directory missing: $RootDir/$Corr"
else
  cd $RootDir/$Corr
  Corr=../$Corr
  # 3pt functions
  for d in 3pt_*
  do
    Spec="${d:4}"
    if [ "${Spec: -2}" = "p2" ]
    then
      Spec="${Spec%p2}"
      pType=p2
      pZero=2_0
      pStar='2_*'
      gammas="gT gXYZ"
    else
      unset pType
      pZero=_0_0_0
      pStar='_*'
      gammas="gT gX gY gZ"
    fi
    Cmd="CRatio --i2 $Corr/2pt${pType}/ --i3 $Corr/$d/ -o $d/"
    Cmd="${Cmd} --efit ,${Spec},g5W0=1000"
    for q in 'h*' l s; do
      Files="'*_${q}_${q}_gT_dt_*_p${pZero}_ps${pZero}_*.h5'"
      echo "$Cmd --type ZV $Files"
    done
    for q in 's_h*' 'l_h*'; do
      for gamma in $gammas; do
        Files="'*_${q}_${gamma}_dt_*_p${pStar}_ps${pStar}.h5'"
        echo "$Cmd --type R $Files"
      done
    done
  done
fi
