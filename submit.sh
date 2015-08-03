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
    rootfilename=$outname"_tree.root"
    rootfilenamedisp="disp_"$outname"_tree.root"
    fileoutname="xmlfiles/scripted_"$outname".xml"
    fileoutnamedisp="xmlfiles/scripted_displaced_"$outname".xml"
#    echo $samplename " " $counter " "$outname" " $infile
    cat dummyconfig.xml | sed -e s%FILENAME%$infile%g -e s%NAME%$outname%g > $fileoutname                                                
    cat dummyconfig.xml | sed -e s%FILENAME%$infile%g -e s%NAME%$outnamedisp%g > $fileoutnamedisp

    echo $fileoutname
    
#    { ./Ntupler $fileoutname >& $sampledir/log_$outname ; mv $rootfilename $sampledir/$rootfilename } &
    { { ./Ntupler $fileoutname ; mv $rootfilename $sampledir/$rootfilename; echo $infile >> done_$sampledir; } & };
    { { ./NtuplerDisplaced $fileoutnamedisp ; mv $rootfilenamedisp $sampledir/$rootfilenamedisp; echo $infile >> done_$sampledir; } & };
    if [ $ii = $njobs ]; then
	echo "now waiting for:" $samplename " " $counter " "$outname" " $infile
	waitabit=$!
	ii=0
	wait $waitabit
    fi
    let "counter++"
    let "ii++"
done
echo "done!" 
