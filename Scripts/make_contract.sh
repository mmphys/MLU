#!/bin/bash
function make_momenta()
{
  momenta=()
  local cond=$2
  local ptot
  let ptot=$1
#  echo All momenta $cond $ptot
  local PX
  local py
  local pz
  let px=-ptot-1
  while (( px++ < ptot )); do
    let py=-ptot-1
    while (( py++ < ptot )); do
      let pz=-ptot-1
      while (( pz++ < ptot )); do
        if (( px * px + py * py + pz * pz $cond ptot )) ; then
#          echo ${px}_${py}_${pz}
          momenta+=(${px}_${py}_${pz})
        fi
      done
    done
  done
  echo ${#momenta[@]} momenta "$cond" $ptot: ${momenta[@]}
}

declare -a momenta
make_momenta 2 ==
make_momenta 2 "<"
make_momenta 1 ==
make_momenta 0 ==
