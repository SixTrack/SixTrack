#!/bin/bash
declare -a error
defaults=$1
add=$2
all=("${defaults[@]}" "${add[@]}")
#simulating map-kinda structure(only bash4.0 natively supports map structure)
declare -i key=0
put()
{
    key=$key+1
    options=("${options[@]}" "$1")
    options_needs=("${options_needs[@]}" "$2")
    options_excludes=("${options_excludes[@]}" "$3")

}

option_check()
{   
    j="$1";
    needed=("${options_needs[${j}]}")
    excluded=("${options_excludes[${j}]}")
    for item in ${needed[@]};do
      if ! echo ${all[@]} | grep -w -q "$item";then
        all=("${all[@]}" "$item")
        add=("${add[@]}" "$item")
      fi
    done

    for item in ${excluded[@]};do
      if echo ${defaults[@]} | grep -w -q "$item";then
        #code to remove conflicting default item from all provided options
        declare -a dummy_all    
        for index in ${all[@]};do
          [[ $index != "$item" ]] && dummy_all+=($index)
        done
        all=("${dummy_all[@]}")
        unset dummy_all
      fi
      if echo ${add[@]} | grep -w -q "$item";then
        error=("error: added option $item (either given by user explicitly or required as dependency) is conflicting with ${options[${i}]}")
      fi 
    done

}


#add new cases below as::  put "option" "needed dependencies" "exclusions"
#put empty parenthesis("") if no exclusions or dependencies(this is important for avoiding worng indexing)

put "hdf5" "collimat" ""
put "hugenblz" "" "bignblz"
put "beamgas" "collimat" "hugenblz bignblz"
put "da" "" "collimat cpss bpm"
put "naglib" "da" ""
put "collimat" "" "da cpss bpm cr"


declare -i i=0
declare -i k=0
while [ $k -lt 10 ];do                         #running test 10 times just to ensure no new conflict arise in items added in 
while [ $i -lt $key ];do                         
 for index in ${all[@]};do
     [[ $index = "${options[$i]}" ]] && option_check "$i"
done
i=$i+1 
done
i=0
k=$k+1
done
if echo ${error} | grep -w -q "error";then
echo ${error}
else
echo ${all[@]}
fi
