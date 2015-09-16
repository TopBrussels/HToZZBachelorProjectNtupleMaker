#!/bin/bash

njobs=3
counter=0
ii=0
samplename=$2
sampledir=$samplename
rm -rf $sampledir
mkdir $sampledir
for infile in `cat $1`
do 
    outname=$samplename"_"$counter
    rootfilenamedisp=$outname"_displacedtree.root"
    xmlfileoutnamedisp="xmlfiles/scripted_displaced_"$outname".xml"
#    echo $samplename " " $counter " "$outname" " $infile
    echo $rootfilenamedisp " " $fileoutnamedisp " " $infile

    cat dummyconfig.xml | sed -e s%FILENAME%$infile%g -e s%OUTNAME%$outname%g > $xmlfileoutnamedisp
#	cat $fileoutnamedisp
#	echo   "./NtuplerDisplaced "$fileoutnamedisp" ; mv "$rootfilenamedisp $sampledir/$rootfilenamedisp
	./NtuplerDisplaced $xmlfileoutnamedisp
	mv $rootfilenamedisp $sampledir/.
#{ { ./NtuplerDisplaced $xmlfileoutnamedisp ; mv $rootfilenamedisp $sampledir/$rootfilenamedisp } & };
#    if [ $ii = $njobs ]; then
#		echo "now waiting for:" $samplename " " $counter " "$outname" " $infile
#		waitabit=$!
#		ii=0
#		wait $waitabit
#    fi
    let "counter++"
    let "ii++"
done
echo "done!" 
