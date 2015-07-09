#!/bin/bash

inputdir="/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v2/"
SKIMNAME=$2
skimfile=$1
TOPPRODPATH=`pwd`/../TopTreeProducer

echo "#usage: . submit_skimmer.sh /path/to/my/skim.xml myskimname"

for dataset in `ls $inputdir`
do
    fullpath=$inputdir$dataset
#    echo $fullpath
    for subdataset in `ls $fullpath`
    do
#        echo $subdataset
	echo " "
        datasetpath=$fullpath"/"$subdataset
	newdatasetpath=$datasetpath/`ls $datasetpath`
	finaldatasetpath=$newdatasetpath/`ls $newdatasetpath`
	datasetpath=$finaldatasetpath
#        echo $datasetpath
        echo "nohup python "`pwd`"/../AutomaticTopTreeProducer/SkimTopTree.py --srmcp -t "$TOPPRODPATH" -l "$datasetpath"/ -s "$skimfile"  --publish "$SKIMNAME"/"$subdataset" -j 4 -n 10 -w 100 --email fblekman@cern.ch  --announce"

#	echo "nohup python /user/fblekman/localgrid/TopTreeWork74X/CMSSW_7_4_3/src/TopBrussels/AutomaticTopTreeProducer/SkimTopTree.py --use-pbs  -t /user/fblekman/localgrid/TopTreeWork74X/CMSSW_7_4_3/src/TopBrussels/TopTreeProducer/ -l "$datasetpath"/ -s /user/fblekman/localgrid/TopTreeWork74X/CMSSW_7_4_3/src/TopBrussels/HToZZBachelorProjectNtupleMaker/skim.xml --publish SkimTests74X/"$subdataset"  -j 4 -n 10 -w 100 --email fblekman@cern.ch  --announce &; sleep 1"

    done
done