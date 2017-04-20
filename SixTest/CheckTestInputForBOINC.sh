#!/usr/bin/env bash
# Script to check that the Sixin.zip files matches the general inputs

set -e

for i in $(ls -d */); do
    if [ ! -f ${i}fort.3 ]; then
	#echo "Folder '"$i"' contains no tests; skipping."
	continue
    fi

    #echo "Checking folder '"${i%%/}"'";
    
    if [ -f ${i}Sixin.zip ]; then
	mkdir ${i}/zip
	cp ${i}/Sixin.zip ${i}/zip
	cd ${i}/zip
	
	unzip -qq Sixin.zip

	for j in $(ls fort.*); do
	    set +e
	    diff -u ${j} ../${j}
	    if [[ $? != 0 ]]; then
		echo "ERROR: Files '"${j}"' does not match in '"${i%%/}"'"
	    fi
	    set -e
	done
	
	if [ -f ../extra_inputs.txt ]; then
	    for k in $(cat ../extra_inputs.txt); do
		set +e
		diff -u ${k} ../${k}
		if [[ $? != 0 ]]; then
		    echo "ERROR: Files '"${j}"' does not match in '"${i%%/}"'"
		fi
		set -e
	    done
	fi
	
	cd ../
	rm -rf zip
	cd ../
    else
	echo "WARNING: No Sixin.zip found in '"${i%%/}"'"
    fi
done
