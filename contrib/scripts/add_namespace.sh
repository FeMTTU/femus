#!/bin/bash

# put this script at the include/ or src/ level (modify .hpp or .cpp accordingly)

# This script is an attempt "to the best you can" of inserting 
# the begin and end of a namespace in the library source files
# The criterion for the BEGIN is: insert it AFTER the LAST #include directive
# The criterion for the END is: insert it BEFORE the LAST #endif directive
# While the first criterion was quite reliable (though still there have been flaws),
# the second criterion basically only works fairly well (but not even exactly) 
# only on .hpp, which have the include guard. For the .cpp it clearly doesn't work at all
# and so the intervention by hand is mandatory.

hand_files=( ) #initialize array of dynamic length
endfile_risky=( )

for file in */*.hpp;
do
sed  -n '/#include/=' $file > tmp_begin.txt;
# cat tmp_begin.txt; 
sed  -n '/#endif/=' $file > tmp_end.txt;
# cat tmp_end.txt; 
echo "=== $file";
export NUM_LINES_TMP_BEGIN=$(wc -l tmp_begin.txt | cut -d " " -f 1); 
export   NUM_LINES_TMP_END=$(wc -l tmp_end.txt   | cut -d " " -f 1); 
# echo "The number of lines with #include directives is " $NUM_LINES_TMP_BEGIN; 
# echo "The number of lines with #endif directives is " $NUM_LINES_TMP_END; 

if (test $NUM_LINES_TMP_BEGIN -gt 0) then #we use the opening criterion to treat files separately

export START_ROW_BEGIN=`tail -n 1 tmp_begin.txt`;
export START_ROW_END=`tail -n 1 tmp_end.txt`;
START_ROW_BEGIN=$(( ${START_ROW_BEGIN} + 2 ))
  START_ROW_END=$(( ${START_ROW_END}   - 2 ))

TOTAL_ROWS=`wc -l $file | cut -d " " -f 1`;
ENDFILE_GAP=$(( $TOTAL_ROWS - $START_ROW_END + 1)); 
echo "The gap between the last endif and the end of file is " $ENDFILE_GAP "*********************"
if (test  $ENDFILE_GAP  -gt 1) 
then
endfile_risky+=($file)
fi
# echo "Start from row " $START_ROW_BEGIN; 
sed  -e  ''$START_ROW_BEGIN'i\\n\nnamespace femus {\n\n'  $file #-i for infile
START_ROW_END=$(( ${START_ROW_END}   + 5 ))  # this is because 5 lines have just been inserted!!!
sed  -e    ''$START_ROW_END'i\\n\n} //end namespace femus\n\n' $file #-i for infile

elif (test $NUM_LINES_TMP_BEGIN -eq 0); 
then echo "*************** $file will be treated by hand ***************" ; 
hand_files+=($file)
fi; 
done

echo "^^^^^^^^^^These are the files to be handled by hand^^^^^^^^^"
for i in $( seq 0 1 ${#hand_files[@]} )
do echo ${hand_files[$i]}
done

echo "^^^^^^^^^^These files might have issues at the end^^^^^^^^^"
for i in $( seq 0 1 ${#endfile_risky[@]} )
do echo ${endfile_risky[$i]}
done


rm tmp_begin.txt
rm tmp_end.txt