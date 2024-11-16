#!/bin/bash

# Default arguments if not supplied
src="$1"
dst="$2"
if [[ "$src" == "" ]]; then src="."; fi
if [[ "$dst" == "" ]]; then dst="smeared"; fi

echo "Moving directories from \"$src\" to \"$dst\", duplicates to \"$dup\""

# Ensure directories exist
if [[ ! -d "$src" ]]
then
  echo "Error: source directory \"$src\" does not exist"
  exit 2
fi
if [[ ! -d "$dst" ]]
then
  echo "Error: destination directory \"$dst\" does not exist"
  exit 2
fi
dup="${dst}_duplicate"
mkdir -p "$dup"

for d in $src/*.3200
do
  if [[ -d "$d" ]]
  then
    this="${d##*/}"
    if [[ -d "$dst/$this" ]]
    then
      cmd="mv $dst/$this $dup/"
      echo "+++ $cmd"
      eval "$cmd"
    fi
    cmd="mv $d $dst/"
    echo "$cmd"
    eval "$cmd"
  fi
done
