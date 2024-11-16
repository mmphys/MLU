#!/usr/bin/env bash

[[ $Random = "" ]] && Random=517671340
[[ $Log = "" ]] && Log=log.txt

function DoFit()
{
  local d=$Base/${1}_
  local Args="$2"
  MultiFit --Hotelling 1e-100 -e 2 -t5:9 $Args -o $d ${Data}_A_{A,B}.fold.$Random.h5 \
    &> ${d}${Base}.$Random.$Log
}

# Sample sizes in synthetic data
for s in 100 1000
do
  # Bootstrap replica counts
  for n in 1000 10000
  do
  if (( s <= n ))
  then
    # Make test data
    Base=S${s}Boot${n}
    mkdir -p $Base
    Data=$Base/$Base
    SynthData -r $Random -s $s -n $n $Data &> $Data.$Log
    
    # Fit test data
    DoFit Corr_Binned "--covsrc binned"
    for CorBoot in 100 1000 10000
    do
      if (( CorBoot <= n ))
      then
        DoFit Corr_Boot$CorBoot "--covsrc bootstrap,$CorBoot"
      fi
    done
  fi
  done
done
