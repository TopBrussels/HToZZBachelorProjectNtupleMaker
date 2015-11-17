rm -rf ntuples/*.root
ls ../ | grep _Ntupler | awk -F"_Ntupler" '{print "hadd ntuples/"$1".root ../"$1"_Ntupler/*.root"}' | sh
