#!/bin/bash

defaults=$1
add=$2
remove=$3
user_given=("${defaults[@]}" "${add[@]}" "${remove[@]}")
list=$4

declare -a unknown

#checking if given options exist in allowed list 
for index in ${user_given[@]};do
  if ! echo ${list[@]} | grep -w -q "$index";then
    unknown+=($index)
  fi
done

if ! [ ${#unknown[@]} -eq 0 ];then
  error=("error - following options dont exist: ${unknown[@]}")
else
  #this means user didnt give unkown option and dependency resolving should be done now
  #simulating map-kinda structure(only bash4.0 natively supports map structure)
  all=("${defaults[@]}" "${add[@]}")
  declare -i key=0
  declare -i key_ifnot=0
  declare -i i=0
  declare -i k=0
  declare -a dummy_all
  declare -a error
  declare -a all_lasttime
  put()
  {
      key=$key+1
      options=("${options[@]}" "$1")
      options_needs=("${options_needs[@]}" "$2")
      options_excludes=("${options_excludes[@]}" "$3")

  }

  put_ifnot()
  {
      key_ifnot=$key_ifnot+1
      options_ifnot=("${options_ifnot[@]}" "$1") ;
      options_ifnot_needs=("${options_ifnot_needs[@]}" "$2")
      options_ifnot_excludes=("${options_ifnot_excludes[@]}" "$3");

  }

  option_check()
  {   
	#this checks if required options in "presence" of certain option are present and excluding option are absent
	#if conflict arises, the option which is in default is removed but if both conflicting option is in add or required by any option in add then gives error 

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

	  for index in ${all[@]};do
	    [[ $index != "$item" ]] && dummy_all+=($index)
	  done
	  all=("${dummy_all[@]}")
	  unset dummy_all
	fi
	if echo ${add[@]} | grep -w -q "$item";then
	  if echo ${defaults[@]} |  grep -w -q "$options[$j]";then
	   for index in ${all[@]};do
	    [[ $index != "$options[$j]" ]] && dummy_all+=($index)
	  done
	  all=("${dummy_all[@]}")
	  unset dummy_all
	  break
	  else
	  error=("error - added option $item (either given by user explicitly or required as dependency) is incompatible with ${options[${i}]}")
	  fi
	fi 
      done

  }

  options_ifnot_check()
  {      # this checks if required options in "absence" of certain option are present and excluding options are absent
	 # if excluding options are in add array, then this gives error

      j="$1";
      needed_ifnot=("${options_ifnot_needs[${j}]}")
      excluded_ifnot=("${options_ifnot_excludes[${j}]}")

      for item in ${needed_ifnot[@]};do
	if ! echo ${all[@]} | grep -w -q "$item";then
	  all=("${all[@]}" "$item")
	  add=("${add[@]}" "$item")
	fi
      done

      for item in ${excluded_ifnot[@]};do
	if echo ${defaults[@]} | grep -w -q "$item";then
	  #code to remove conflicting default item from all provided options
	  for index in ${all[@]};do
	    [[ $index != "$item" ]] && dummy_all+=($index)
	  done
	  all=("${dummy_all[@]}")
	  unset dummy_all
	fi
	if echo ${add[@]} | grep -w -q "$item";then
	  error=("error - added option $item (either given by user explicitly or required as dependency) is incompatible with absence of ${options_ifnot[${i}]}")
	fi 
      done

  }

  ## LIST OF DEPENDENCIES: ##
  
  # Add new cases below as::  put "option" "needed dependencies" "exclusions"
  # put empty parenthesis("") if no exclusions or dependencies(this is important for avoiding worng indexing)
  put "hdf5" "collimat" ""
  put "bonic" "cpss crlibm" ""
  put "api" "boinc" ""
  put "beamgas" "collimat" "bignblz hugenblz"
  put "hugenblz" "" "bignblz"
  put "da" "" "collimat cpss bpm"
  put "collimat" "" "da cpss bpm cr crlibm"
  put "cpss" "crlibm cr" "cernlib" 
  put "m64" "" "ifort nagfor pgf90 g95 lf95 cernlib bonic m32 naglib"
  put "beamgas" "collimat" ""
  
  # Add below the case which are to be followed in absence of certain options
  # as put_ifnot "option" "needed dependencies" "exclusions"
  # i.e. 'put_ifnot "da" "" "naglib"' means "if not da, then naglib is forbidden",
  # and 'put_ifnot "m64" "m32" ""' means "if not m64, then m32 is required".
  put_ifnot "da" "" "naglib"
  put_ifnot "m64" "m32" ""

  while [[ $k -lt 10  && ${all_lasttime[@]} != ${all[@]} ]];do             #running test max 10 times just to ensure no new conflict arise in items added
  unset error
  all_lasttime=${all[@]}
  while [ $i -lt $key ];do                         
   for index in ${all[@]};do
       [[ $index = "${options[$i]}" ]] && option_check "$i"
  done
  i=$i+1 
  done
  i=0
  while [ $i -lt $key_ifnot ];do                         
   if ! echo ${all[@]} | grep -w -q "${options_ifnot[$i]}";then
    options_ifnot_check "$i"
   fi
  i=$i+1 
  done
  i=0
  k=$k+1
  done
fi

#finally providing make with either error or correct options
if ! [ ${#error[@]} -eq 0 ];then
  echo ${error[@]}
else
  echo ${all[@]}
fi
