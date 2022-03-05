#!/bin/bash

set -e

jupyter nbconvert --to script 'generate_configuration_gromacs.ipynb'

for cions in 1 2 3 4 5;
do
	newline='c = '$cions
	linetoreplace=$(cat generate_configuration_gromacs.py | grep 'c =')
	sed -i '/'"$linetoreplace"'/c\'"$newline" generate_configuration_gromacs.py

 	DIRNAME=cs_${cions}M
	if [ ! -d "../$DIRNAME" ]; then
		echo "creating folder"
		mkdir ../$DIRNAME	
	fi

	python3 generate_configuration_gromacs.py

	mv conf.gro topol.top ../$DIRNAME
done

rm generate_configuration_gromacs.py
