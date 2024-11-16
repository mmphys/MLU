#!/bin/bash
case $HOSTNAME in
  "S44-40EGHV2F") src="$HOME/src/Utility/run";;
  *)              src="$HOME/../dc-erbe1/contract/nvec_studies";;
esac

# perform a bootstrap on data in $src/$1, optionally excluding $2
function do_bootstrap {
  cmd="bootstrap -o $1.@corr@.@type@ $src/$1/*.h5"
  if [ ! -z "$2" ] ; then cmd="$cmd -i $src/$1/$2"; fi
  #echo $cmd
  eval $cmd
}

for batch in light_100 light_75 light_50; do
  cmd="do_bootstrap $batch"
  if [ "$batch" = "light_50" ] ; then cmd="$cmd phiL_0_rhoL.3000.h5"; fi
  eval $cmd
done
#bootstrap -o light_100.@corr@.@type@ $src/light_100/*.h5
#bootstrap -o light_75.@corr@.@type@  $src/light_75/*.h5
#bootstrap -o light_50.@corr@.@type@  $src/light_50/*.h5 -i $src/light_50/phiL_0_rhoL.3000.h5
