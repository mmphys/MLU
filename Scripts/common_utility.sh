#!/usr/bin/env bash

PlotPathSplit()
{
  unset mmplotfile_name_no_ext
  unset mmplotfile_base
  unset mmplotfile_split
  unset mmplotfile_type
  unset mmplotfile_seed
  unset mmplotfile_ext
  unset mmplotfile_corr
  unset mmplotfile_corr_all
  unset mmplotfile_ti
  unset mmplotfile_tf
  unset mmplotfile_ops
  unset mmplotfile_ops_all
  # Split into path and name
  mmplotfile_name="$1"
  mmplotfile_path="${mmplotfile_name%/*}"
  if [[ "$mmplotfile_path" == "$mmplotfile_name" ]]
  then
    mmplotfile_path=
  else
    mmplotfile_path="$mmplotfile_path/"
    mmplotfile_name="${mmplotfile_name##*/}"
  fi
  # Break filename into components
  mmplotfile_split=(${mmplotfile_name//./ })
  # Save specific components
  if (( ${#mmplotfile_split[@]} < 2 ))
  then
    mmplotfile_name_no_ext="${mmplotfile_name}"
    mmplotfile_base="${mmplotfile_name}"
  else
    mmplotfile_name_no_ext="${mmplotfile_name%.*}"
    mmplotfile_ext="${mmplotfile_split[@]: -1:1}"
    if (( ${#mmplotfile_split[@]} < 4 ))
    then
      mmplotfile_base="${mmplotfile_split[@]:0:$(( ${#mmplotfile_split[@]} - 1 ))}"
    else
      mmplotfile_type="${mmplotfile_split[@]: -3:1}"
      mmplotfile_seed="${mmplotfile_split[@]: -2:1}"
      # If it's six words long, then dig out extra info
      if (( ${2:-4} < 6 || ${#mmplotfile_split[@]} < 6 ))
      then
        mmplotfile_base="${mmplotfile_split[@]:0:$(( ${#mmplotfile_split[@]} - 3 ))}"
      else
        mmplotfile_base="${mmplotfile_split[@]:0:$(( ${#mmplotfile_split[@]} - 5 ))}"
        mmplotfile_ops_all="${mmplotfile_split[@]: -4:1}"
        mmplotfile_ops=(${mmplotfile_ops_all//_/ })
        mmplotfile_corr_all="${mmplotfile_split[@]: -5:1}"
        local mmplotfile_tmp=(${mmplotfile_corr_all//_/ })
        if (( ${#mmplotfile_tmp[@]} != 3 ))
        then
          mmplotfile_corr="${mmplotfile_corr_all}"
        else
          mmplotfile_corr="${mmplotfile_tmp[@]: 0:1}"
          mmplotfile_ti="${mmplotfile_tmp[@]: 1:1}"
          mmplotfile_tf="${mmplotfile_tmp[@]: 2:1}"
        fi
      fi
      mmplotfile_base="${mmplotfile_base// /.}"
    fi
  fi
}

export -f PlotPathSplit
export mmplotvar_txt=txt # default extension for data files in text
