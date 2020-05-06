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
  echo "wall Wall time (in hours) for each file (default: 48)"
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
let wall_hours=${wall:=48}
if (( $? || wall_hours <= 0 || wall_hours > 48)); then
  echo "Wall hours must be between 1 and 48"
  exit
fi

#echo "start=$start, NumFiles=$NumFiles, Inc=$Inc"
#echo "num=$num, step=$step"

bOK=1
template="template"
base="${1%.}" # Get rid of a trailing full-stop if present
if   [[ -e $base.$template.sh  ]]; then jobext="sh" ; jobDescr="GPU"
elif [[ -e $base.$template.pbs ]]; then jobext="pbs"; jobDescr="CPU"
else
  bOK=0
  echo "Neither $base.$template.sh nor $base.$template.pbs found"
fi
if [[ ! -e $base.$template.xml ]]; then bOK=0; echo "$base.$template.xml missing"; fi
if (( ! bOK )); then exit; fi

templscript=$base.$template.$jobext
templxml=$base.$template.xml
echo "Making $NumFiles $jobDescr jobs"

for (( ; NumFiles-- ; start += Inc )); do
  end=$(( start + ( num - 1 ) * step + 1 ))
  jobscript=$base.$start.$jobext
  jobxml=$base.$start.xml
  sed -e "s|@start@|$start|g" -e "s|@end@|$end|g" -e "s|@step@|$step|g" -e "s|@xml@|$jobxml|g" -e "s|@wall@|$wall|g" $templscript > $jobscript
  sed -e "s|@start@|$start|g" -e "s|@end@|$end|g" -e "s|@step@|$step|g" $templxml > $jobxml
done
