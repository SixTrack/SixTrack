#!/bin/bash
for i in $(ls -d */); do
#echo ${i%%/};
if [ -f ${i}Sixin.zip ]; then
	mkdir ${i}/zip
	cp ${i}/Sixin.zip ${i}/zip
	cd ${i}/zip
	unzip -qq Sixin.zip
	for j in $(ls fort.*); do
		diff -u ${j} ../${j}
	done
	cd ../
	rm -rf zip
	cd ../
fi
done
