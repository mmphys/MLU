#!/bin/bash

if (( $# < 2 || $# > 4 )); then
  echo "Create n copies of .pbs and corresponding .xml from template"
  echo "\$1 Name of template files (without .pbs/.xml)"
  echo "\$2 Number of first configuration"
  echo "\$3 Number of job files to create (default: 1)"
  echo "\$4 Increment between job files   (default: 40)"
  echo "Optional environment variables to set prior to calling:"
  echo "num  Number of configurations each file should process (default: 1)"
  echo "step Step between configurations in each file (default: 40)"
  exit
fi

num=${num:-1}
step=${step:-40}

let start=$2
if (( $? || start <= 0)); then echo "First configuration must be a number > 0"; exit; fi
let NumFiles=${3:-1}
if (( $? || NumFiles <= 0)); then echo "Number of job files must be a number > 0"; exit; fi
let Inc=${4:-40}
if (( $? || Inc <= 0)); then echo "Increment must be a number > 0"; exit; fi

echo "start=$start, NumFiles=$NumFiles, Inc=$Inc"
echo "num=$num, step=$step"

bOK=1
template="${1%.}"
if [[ ! -e $template.pbs ]]; then bOK=0; echo "$template.pbs missing"; fi
if [[ ! -e $template.xml ]]; then bOK=0; echo "$template.xml missing"; fi
if (( ! bOK )); then exit; fi

export start
export step

for (( ; NumFiles-- ; start += Inc )); do
  export end=$(( start + ( num - 1 ) * step + 1 ))
  base=${template//.*}.$start
  export xml=$base.xml
  #envsubst < $template.pbs > $base.pbs
  #envsubst < $template.xml > $xml
  sed -e "s|@start@|$start|g" -e "s|@end@|$end|g" -e "s|@step@|$step|g" -e "s|@xml@|$xml|g" $template.pbs > $base.pbs
  sed -e "s|@start@|$start|g" -e "s|@end@|$end|g" -e "s|@step@|$step|g" $template.xml > $xml
done
