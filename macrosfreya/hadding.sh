#rm -rf ntuples/*.root
ls ../ | grep _Ntupler | awk -F"_Ntupler" '{print "rm -rf ntuples/"$1".root; hadd ntuples/"$1".root ../"$1"_Ntupler/*.root"}' | grep Run2015 | sh
ls ../ | grep _Ntupler | awk -F"_Ntupler" '{print "rm -rf ntuples/"$1".root; hadd ntuples/"$1".root ../"$1"_Ntupler/*.root &"}' | grep 13TeV | sh
