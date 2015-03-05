#!/bin/bash

# Author: G. Bornia

# This script updates all include guards to be consistent with the filename (including directory structure),
# so as to avoid any possibility of equal include guards that would result
# This script should be run every time a new hpp file is added or you change a directory name

suspect_threshold=30   #if there is no #ifndef before this line number,we suspect that the include guard was forgotten

folder_under_exam=`readlink -f $1`

    if  [ "${folder_under_exam}" = "" ]; then
      echo "$(tput setaf 1)Specify the folder inside which you want to update all include guards in the header files$(tput sgr 0)"
      exit
    fi


count=0
for file in $(find ${folder_under_exam} -iname *.hpp)
do linenum_ifndef=`grep -n -m 1 "#ifndef"  $file | cut -d ":" -f1`
  echo The line with the beginning of the include guard to be updated is ${linenum_ifndef}
  
    if  [ "${linenum_ifndef}" = "" ]; then
      echo The header file ${file} surely misses its include guard
      exit
    fi
    
    if   [ ${linenum_ifndef} -gt ${suspect_threshold} ]; then
      echo The header file ${file} might miss its include guard
      exit
    fi
  
  closest_folder_holding_file=$(basename $(dirname $file))
  file_name_simple=$(basename $file .hpp)
  guard_var=__femus_${closest_folder_holding_file}_${file_name_simple}_hpp__
   first_line=$(echo "#ifndef" ${guard_var})
  second_line=$(echo "#define" ${guard_var})
  echo ${first_line}
  echo ${second_line}
    sed "${linenum_ifndef} c\\${first_line}" -i $file
    sed "$((linenum_ifndef+1)) c\\${second_line}" -i $file
  count=$((count + 1))
done
 
echo Number of header files in the folder ${folder_under_exam}: $count