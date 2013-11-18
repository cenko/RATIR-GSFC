#!/bin/bash

rm *.list

for dir in * ; do 
    if [ -d "$dir" ]; then
	file0='C0_'${dir#./}'.list'
	file1='C1_'${dir#./}'.list'
	file2='C2_'${dir#./}'.list'
	file3='C3_'${dir#./}'.list'
	cd $dir
	for f in *C0*.fits ; do
	    echo $f >> '../'$file0
	done
	for f in *C1*.fits ; do
	    echo $f >> '../'$file1
	done
	for f in *C2*.fits ; do
	    echo $f >> '../'$file2
	done
	for f in *C3*.fits ; do
	    echo $f >> '../'$file3
	done
	cd ../
    fi
done

cd ../
