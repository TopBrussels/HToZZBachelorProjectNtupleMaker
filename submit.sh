#!/bin/bash

nsplit=10
njobs=25
counter=0
ii=0
samplename=$2
sampledir=$samplename"_Ntupler"
rm -rf $sampledir
mkdir $sampledir

samplefile=$sampledir/workfile.txt
rm -rf $samplefile
jj=0
outname=""
for infile in `cat $1`
do
    if [[ $jj != 0 ]]
    then
	outname=$outname,
    fi
    echo "now : "$jj " files, next: "$infile
    outname=$outname'dcap://maite.iihe.ac.be'$infile
#    echo $outname
    let "jj++"
    if [[ $jj == $nsplit ]]
    then
        jj=0
	echo $outname >> $samplefile
	outname=""
    fi
done
echo " Samples: "
echo $outname >> $samplefile
cat $samplefile
echo "--------------"

for infile in `cat $samplefile`
do 
    outname=$samplename"_"$counter
    rootfilename=$outname"_tree.root"
    rootfilenamedisp="disp_"$outname"_tree.root"
    fileoutname="xmlfiles/scripted_"$outname".xml"
    fileoutnamedisp="xmlfiles/scripted_displaced_"$outname".xml"
#    echo $samplename " " $counter " "$outname" " $infile
    cat dummyconfig.xml | sed -e s%FILENAME%$infile%g -e s%NAME%$outname%g > $fileoutname                                                

#    echo $fileoutname

    pbsoutname=$sampledir"/"$outname".pbs"
    echo "pbs file = "$pbsoutname
    cat dummypbs.pbs | sed -e s%working%`pwd`%g > $pbsoutname
    echo "\$DESTDIR/Ntupler \$DESTDIR/"$fileoutname >> $pbsoutname
    echo " mv "$rootfilename" \$DESTDIR/"$sampledir"/." >> $pbsoutname
    queue="localgrid"
    cat $pbsoutname | sed -e s%"walltime=1:14:59"%"walltime=0:30:00"%g > tmpfile
    mv tmpfile $pbsoutname
    if [[ `echo $outname | grep 13TeV` ]]  # it's MC!
    then
	queue="localgrid"
	# and increase the run time:
	cat $pbsoutname | sed -e s%"walltime=0:30:00"%"walltime=4:00:00"%g > tmpfile
    mv tmpfile $pbsoutname
    fi


    qsub -q $queue $pbsoutname

    if [ "$?" = "1" ]; then
	#failed to submit!
	echo "failed to submit job "$pbsoutname", trying again in 10 sec"
	sleep 10
	qsub -q $queue $pbsoutname
	if [ "$?" = "1" ]; then
	    sleep 1
            qsub -q $queue $pbsoutname
	fi
    fi


    let "counter++"
    let "ii++"
    if [[ $ii == $njobs ]]
    then
	echo "sleeping for a while to check batch privileges"
	sleep 30
	ii=0
    fi
done
echo "done!" 
