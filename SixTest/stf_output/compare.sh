#!/bin/bash
cd stf/
gfortran new_converter.f
./a.out
cd ../not_stf
gfortran old_converter.f
./a.out
cd ../ 
diff stf/output.dat not_stf/older.dat
var=$?
if test $var -eq 0
then
echo "output is SAME"
else
echo "output is  DIFFERENT"
fi
cd stf; rm output.dat a.out
cd ../not_stf/;rm older.dat a.out
cd ../
