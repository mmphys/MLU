#!/bin/bash

# Copy Grid correlators from source to destination, dropping all other metadata
dest=$1
source=${@:2}
#save pattern matching state before turning on extended matches
extglobstate="$(shopt -p extglob)"
shopt -s extglob
# Get rid of consecutive / in destination
pattern="+(/)"
string=/
dest="${dest//$pattern/$string}"
# Get rid of trailing / in destination
if [[ ${#dest} > 1 ]]; then
  dest="${dest%/}"
fi

# Validate arguments
if [[ "$source" == "" ]]; then
  echo "Second and subsequent arguments should be the source(s), e.g. *.h5"
elif [[ "$dest" == "" ]]; then
  echo "First argument should be the destination directory"
elif [ ! -w "$dest" ]; then
  echo "$dest doesn't exist/not writeable"
else
  echo copying $source to $dest
  for file in $source; do
    base=${file%.*}
    base=${base##*/}
    group=${base%.*}
    cmd="h5copy -i $file -o $dest/$base.tiny.h5 -s $group/correlator -d correlator"
    echo $cmd
    eval $cmd
  done
fi
eval $extglobstate
