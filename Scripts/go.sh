#!/bin/bash

if (( $# < 2 || $# > 4 )); then
  echo "Create n copies of .pbs and corresponding .xml from template"
  echo "\$1 Name of template file"
  echo "\$2 Number of first configuration"
  echo "\$3 Number of job files to create (default: 1)"
  echo "\$4 Increment between job files   (default: 40)"
  echo "Optional environment variables to set prior to calling:"
  echo "num  Number of configurations each file should process (default: 1)"
  echo "step Step between configurations in each file (default: 40)"
  echo "wall Wall time (in hours) for each file (default: 48/24 CPU/GPU)"
  exit
fi

template="$1"
jobext=${template##*.}
base=${template%.*}
templateXml=${base}.xml
base=${base%.*}
shopt -s nocasematch
case $jobext in
  sh) jobDescr="GPU"; walldefault=24;;
  pbs) jobDescr="CPU"; walldefault=48;;
  *) echo "extension should be .pbs for CPU or .sh for GPU"; exit 1;;
esac
bOK=1
if [[ ! -e $template    ]]; then bOK=0; echo "$template missing"; fi
if [[ ! -e $templateXml ]]; then bOK=0; echo "$templateXml missing"; fi
if (( ! bOK )); then exit; fi

let start=$2
if (( $? || start <= 0)); then echo "First configuration must be a number > 0"; exit; fi
let NumFiles=${3:-1}
if (( $? || NumFiles <= 0)); then echo "Number of job files must be a number > 0"; exit; fi
let Inc=${4:-40}
if (( $? || Inc <= 0)); then echo "Increment must be a number > 0"; exit; fi
let wall_hours=${wall:-${walldefault}}
if (( $? || wall_hours <= 0 || wall_hours > 48)); then
  echo "Wall hours must be between 1 and 48"
  exit
fi

num=${num:-1}
step=${step:-40}

#echo "start=$start, NumFiles=$NumFiles, Inc=$Inc"
#echo "num=$num, step=$step"

echo "Making $NumFiles $jobDescr jobs"

for (( ; NumFiles-- ; start += Inc )); do
  end=$(( start + ( num - 1 ) * step + 1 ))
  jobscript=$base.$start.$jobext
  jobxml=$base.$start.xml
  #Echo making $jobscript and $jobxml
  sed -e "s|@start@|$start|g" -e "s|@end@|$end|g" -e "s|@step@|$step|g" -e "s|@xml@|$jobxml|g" -e "s|@wall@|${wall_hours}|g" $template > $jobscript
  sed -e "s|@start@|$start|g" -e "s|@end@|$end|g" -e "s|@step@|$step|g" $templateXml > $jobxml
done
