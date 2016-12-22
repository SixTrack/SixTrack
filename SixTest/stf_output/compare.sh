#!/bin/bash
# Author: Vikas Gupta , email-id: vikasgupta.rmps@gmail.com
# script to read in output from stf and non-stf sixtrack run and compare them 
# for any differences 
not_stf_path=$1
stf_path=$2
cd stf; rm -f output.dat new_converter.out fort.*
cd ../not_stf/;rm -f older.dat old_converter.out fort.*
cd ../
cd not_stf
cp -r $not_stf_path/fort.* .
gfortran old_converter.f -o old_converter.out -Wall
./old_converter.out
cd ../stf/
cp $stf_path .
gfortran new_converter.f -o new_converter.out -Wall
./new_converter.out
cd ../
diff stf/output.dat not_stf/older.dat
var=$?
if test $var -eq 0
then
echo "output is SAME"
else
echo "output is  DIFFERENT"
fi
cd stf; rm -f output.dat new_converter.out fort.*
cd ../not_stf/;rm -f older.dat old_converter.out fort.*
cd ../
