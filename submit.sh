#!/bin/bash

njobs=3
counter=0
ii=0
samplename=$2
sampledir=$samplename"_Ntupler"
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

    echo $fileoutname

    pbsoutname=$sampledir"/"$outname".pbs"
    echo "pbs file = "$pbsoutname
    cat dummypbs.pbs > $pbsoutname
#	echo "cp \$DESTDIR/"$fileoutname" infile.xml" >> $pbsoutname
    echo "\$DESTDIR/Ntupler \$DESTDIR/"$fileoutname >> $pbsoutname
    echo " mv "$rootfilename" \$DESTDIR/"$sampledir"/." >> $pbsoutname
    cat $pbsoutname | sed -e s%"/user/fblekman/localgrid/"%"/localgrid/fblekman/"%g > tmpfile
    mv tmpfile $pbsoutname
#cat $pbsoutname
    queue="short"
    if [[ `echo $outname | grep 13TeV` ]]  # it's MC!
    then
	queue="localgrid"
	# and increase the run time:
	cat $pbsoutname | sed -e s%"walltime=1:14:59"%"walltime=3:00:00"%g > tmpfile
    mv tmpfile $pbsoutname
    fi

    qsub -q $queue $pbsoutname


    let "counter++"
    let "ii++"
    if [[ $ii == 200 ]]
    then
	echo "sleeping for a while to check batch privileges"
	sleep 30
	$ii=0
    fi
done
echo "done!" 
