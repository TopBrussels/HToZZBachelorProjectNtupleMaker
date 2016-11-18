python skimmer.py -i "../../CMSSW_80X_v1-ntuples/crab_DoubleMuon-Run2016*.root" -o SkimmedDoubleMuon.root -nMu 2 -pt 30 &
python skimmer.py -i "../../CMSSW_80X_v1-ntuples/crab_DoubleEG-Run2016*.root" -o SkimmedDoubleElectron.root -nEl 2 -pt 30 &
python skimmer.py -i "../../CMSSW_80X_v1-ntuples/DYJetsToLL_M-50*.root" -o SkimmedDoubleMuon_DYJetsToLL_M-50.root -nMu 2 -pt 30 -x 5765.4 &
python skimmer.py -i "../../CMSSW_80X_v1-ntuples/DYJetsToLL_M-50*.root" -o SkimmedDoubleElectron_DYJetsToLL_M-50.root -nEl 2 -pt 30 -x 5765.4 &
python skimmer.py -i "../../TT_TuneCUETP8M1_13TeV-powheg-pythia8_*.root" -o SkimmedDoubleElectron_TTpowheg.root -nEl 2 -pt 30 -x 831 &
python skimmer.py -i "../../TT_TuneCUETP8M1_13TeV-powheg-pythia8_*.root" -o SkimmedDoubleMuon_TTpowheg.root -nMu 2 -pt 30 -x 831 & 
