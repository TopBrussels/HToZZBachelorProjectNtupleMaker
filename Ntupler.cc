///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Ntuple.cc: This macro is intended to be an example analysis macro which works out of the box.           /////
/////       It should serve as the first port of call for new users of the TopBrussels framework.             /////
/////      (in addition it is used by Freya for occasional studies when she has time)                         /////
/////     Last Modified: Mon 16 February 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{
    
    clock_t start = clock();
    
    
    /////////////////////
    // Configuration
    /////////////////////
    
    //xml file
    string xmlFileName ="myhiggsconfig.xml";
    
    if (argc > 1)
        xmlFileName = (string)argv[1];
    
    const char *xmlfile = xmlFileName.c_str();
    
    cout << "********************************************************" << endl;
    cout << "used config file: " << xmlfile << endl;
    cout << "********************************************************" << endl;
    
    
    
    cout << "********************************************************" << endl;
    cout<<"creating datasets ..."<<endl;
    cout << "********************************************************" << endl;
    
    //Configuration output format
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    
    ////////////////////////////////////
    /// AnalysisEnvironment
    ////////////////////////////////////
    
    AnalysisEnvironment anaEnv;
    cout << "********************************************************" << endl;
    cout<<"Loading environment ..."<<endl;
    cout << "********************************************************" << endl;
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    
    cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
    
    int verbose = 2;//anaEnv.Verbose;
    
    
    cout << "now done creating AnalysisEnvironmentLoader" << endl;
    cout << "********************************************************" << endl;
    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
    verbose = anaEnv.Verbose;
    float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
    
    Double_t muoneffaverage[2]={0,0};
    Double_t jeteffaverage[2]={0,0};
    Double_t eleeffaverage[2]={0,0};
    
    TLorentzVector worker;
    TTreeLoader treeLoader;
    //    cout << " - Load datasets ..." << endl;
    vector < Dataset* > datasets;
    
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout << "now loaded " << datasets.size() << " datasets" << endl;
    for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
    
    float Luminosity = oldLuminosity;
    string dataSetName="unknown";
    cout << "********************************************************" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
        
        dataSetName = datasets[d]->Name();
        cout << "datasets: " << dataSetName << endl;
    }
    cout << "********************************************************" << endl;
    
    
    
    //Global variable
    //TRootEvent* event = 0;
    
    //nof selected events
    double NEvtsData = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
    ////////////////////////////////////
    //	Loop on datasets
    ////////////////////////////////////
    
    // trigger infos: (from https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger )
    
    
    
    
    std::vector<std::string> mujetstriggers;
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_v2");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v2");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_v2");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v2");
    std::vector<std::string> ejetstriggers;
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v1");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v1");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v1");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v1");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1");
    
    std::vector<std::string> eetriggers;
    eetriggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2");
    eetriggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1");
    
    std::vector<std::string> mumutriggers;
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2");
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2");
    mumutriggers.push_back("HLT_IsoMu20_eta2p1_v2");
    mumutriggers.push_back("HLT_IsoMu20_eta2p1_v1");
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1");
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1");
    
    std::vector<std::string> emutriggers;
    emutriggers.push_back("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2");
    emutriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2");
    emutriggers.push_back("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1");
    emutriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1");
    
    std::vector<std::string> displtriggers;
    displtriggers.push_back("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL_v2");
    displtriggers.push_back("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL_v2");
    
    
    std::map<std::string,std::pair<int,bool> > triggermap;
    // book all these in the trigger map so it can be used later:
    for(UInt_t itrig=0; itrig<mumutriggers.size(); itrig++){
        triggermap[mumutriggers[itrig]]=std::pair<int,bool>(-999,false);
    }
    for(UInt_t itrig=0; itrig<mujetstriggers.size(); itrig++){
        triggermap[mujetstriggers[itrig]]=std::pair<int,bool>(-999,false);
    }
    for(UInt_t itrig=0; itrig<emutriggers.size(); itrig++){
        triggermap[emutriggers[itrig]]=std::pair<int,bool>(-999,false);
    }
    for(UInt_t itrig=0; itrig<eetriggers.size(); itrig++){
        triggermap[eetriggers[itrig]]=std::pair<int,bool>(-999,false);
    }
    for(UInt_t itrig=0; itrig<ejetstriggers.size(); itrig++){
        triggermap[ejetstriggers[itrig]]=std::pair<int,bool>(-999,false);
    }
    for(UInt_t itrig=0; itrig<displtriggers.size(); itrig++){
        triggermap[displtriggers[itrig]]=std::pair<int,bool>(-999,false);
    }
    
    for(std::map<std::string,std::pair<int,bool> >::iterator trigiter = triggermap.begin(); trigiter != triggermap.end(); trigiter++){
        std::pair<int,bool> bla = trigiter->second;
        std::string bla2 = trigiter->first;
        
    }
    
    
    
    
    LumiReWeighting LumiWeights;
    LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt50ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_50ns-Run246908-251883Cert/nominal.root", "pileup", "pileup");
    
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        string previousFilename = "";
        int iFile = -1;
        string dataSetName = datasets[d]->Name();
        
        Int_t isdata;
        
        cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
        if (verbose > 1)
            std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
        //         make root tree file name
        string roottreename = datasets[d]->Name();
        roottreename+="_tree.root";
        //        cout << "creating tree in file " << roottreename << endl;
        
        // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
        fileout->cd();
        
        //////////////////////////////
        // My tree - variables      //
        //////////////////////////////
        Int_t nElectrons;
        Double_t pT_electron[10];
        Double_t phi_electron[10];
        Double_t eta_electron[10];
        Double_t E_electron[10];
        Double_t d0_electron[10];
        Double_t pfIso_electron[10];
        Int_t charge_electron[10];
        Int_t loose_electron[10];
        Int_t medium_electron[10];
        Int_t tight_electron[10];
        Int_t mediumMVA_electron[10];
        Int_t tightMVA_electron[10];
        
        Int_t nMuons;
        Double_t pT_muon[10];
        Double_t phi_muon[10];
        Double_t eta_muon[10];
        Double_t E_muon[10];
        Double_t d0_muon[10];
        Double_t pfIso_muon[10];
        Int_t charge_muon[10];
        
        Int_t nJets;
        Double_t pT_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Double_t dx_jet[20];
        Double_t dy_jet[20];
        Double_t btag_jet[20];
        Double_t missingEt;
        // various weights
        Double_t pu_weight;
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
        Int_t trig_dilepton_emu;
        Int_t trig_dilepton_mumu;
        Int_t trig_dilepton_ee;
        Int_t trig_eplusjets;
        Int_t trig_muplusjets;
        Int_t trig_displaced;
        
        TTree *bookkeeping = new TTree("startevents","startevents");
        bookkeeping->Branch("run_num",&run_num,"run_num/I");
        bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
        bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
        bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
        bookkeeping->Branch("npu",&npu,"npu/I");
        
        // define the output tree
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isdata",&isdata,"isdata/I");
        myTree->Branch("run_num",&run_num,"run_num/I");
        myTree->Branch("evt_num",&evt_num,"evt_num/I");
        myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        myTree->Branch("nvtx",&nvtx,"nvtx/I");
        myTree->Branch("npu",&npu,"npu/I");
        myTree->Branch("trig_dilepton_emu",&trig_dilepton_emu,"trig_dilepton_emu/I");
        myTree->Branch("trig_dilepton_ee",&trig_dilepton_ee,"trig_dilepton_ee/I");
        myTree->Branch("trig_dilepton_mumu",&trig_dilepton_mumu,"trig_dilepton_mumu/I");
        myTree->Branch("trig_eplusjets",&trig_eplusjets,"trig_eplusjets/I");
        myTree->Branch("trig_muplusjets",&trig_muplusjets,"trig_muplusjets/I");
        myTree->Branch("trig_displaced",&trig_displaced,"trig_displaced/I");
        
        
        
        
        
        
        
        
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        myTree->Branch("pT_electron",pT_electron,"pT_electron[nElectrons]/D");
        myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        myTree->Branch("loose_electron",loose_electron,"loose_electron[nElectrons]/I");
        myTree->Branch("medium_electron",medium_electron,"medium_electron[nElectrons]/I");
        myTree->Branch("tight_electron",tight_electron,"tight_electron[nElectrons]/I");
        //        myTree->Branch("mediumMVA_electron",mediumMVA_electron,"mediumMVA_electron[nElectrons]/I");
        //        myTree->Branch("tightMVA_electron",tightMVA_electron,"tightMVA_electron[nElectrons]/I");
        //
        
        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pT_muon",pT_muon,"pT_muon[nMuons]/D");
        myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        
        myTree->Branch("nJets",&nJets, "nJets/I");
        myTree->Branch("pT_jet",pT_jet,"pT_jet[nJets]/D");
        myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        myTree->Branch("dx_jet",dx_jet,"dx_jet[nJets]/D");
        myTree->Branch("dy_jet",dy_jet,"dy_jet[nJets]/D");
        myTree->Branch("btag_jet",btag_jet,"btag_jet[nJets]/D");
        
        myTree->Branch("missingEt",&missingEt,"missingEt/D");
        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
        
        TTree *noselTree = new TTree("noselTree","noselTree");
        noselTree->Branch("isdata",&isdata,"isdata/I");
        noselTree->Branch("run_num",&run_num,"run_num/I");
        noselTree->Branch("evt_num",&evt_num,"evt_num/I");
        noselTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        noselTree->Branch("nvtx",&nvtx,"nvtx/I");
        
        
        noselTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        noselTree->Branch("pT_electron",pT_electron,"pT_electron[nElectrons]/D");
        noselTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        noselTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        noselTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        noselTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        noselTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        noselTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        
        
        noselTree->Branch("nMuons",&nMuons, "nMuons/I");
        noselTree->Branch("pT_muon",pT_muon,"pT_muon[nMuons]/D");
        noselTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        noselTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        noselTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        noselTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        noselTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        noselTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        
        noselTree->Branch("nJets",&nJets, "nJets/I");
        noselTree->Branch("pT_jet",pT_jet,"pT_jet[nJets]/D");
        noselTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        noselTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        noselTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        noselTree->Branch("dx_jet",dx_jet,"dx_jet[nJets]/D");
        noselTree->Branch("dy_jet",dy_jet,"dy_jet[nJets]/D");
        noselTree->Branch("btag_jet",btag_jet,"btag_jet[nJets]/D");
        
        
        noselTree->Branch("missingEt",&missingEt,"missingEt/D");
        noselTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
        
        
        //        myTree->Print();
        
        
        
        
        //open files and load
        treeLoader.LoadDataset (datasets[d], anaEnv);
        
        
        //////////////////////////////////////////////////
        /// Initialize Jet energy correction factors   ///
        //////////////////////////////////////////////////
        
        vector<JetCorrectorParameters> vCorrParam;
        
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
        
        // true means redo also the L1 corrections (see CMS documentation to learn what this means)
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
        
        
        ////////////////////////////////////
        //	Loop on events
        ////////////////////////////////////
        
        // some bookkeeping variables
        nEvents[d] = 0;
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
        
        // some printout
        cout << "running over " << datasets[d]->NofEvtsToRunOver() << endl;
        
        // start event loop
        for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++) // event loop
            //for (unsigned int ievt = 0; ievt < 20000; ievt++) // if lazy for testing
        {
            
            // the objects loaded in each event
            vector < TRootHLTInfo> hltinfo;
            vector < TRootVertex* > vertex;
            vector < TRootMuon* > init_muons;
            vector < TRootElectron* > init_electrons;
            vector < TRootJet* > init_jets_corrected;
            vector < TRootJet* > init_jets;
            vector < TRootMET* > mets;
            vector < TRootGenJet* > genjets;
            
            nEvents[d]++;
            
            if(ievt%1000 == 0)
                std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected" << flush<<"\r";
            
            ////////////////
            // LOAD EVENT //
            ////////////////
            
            TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
            std::map<std::string, std::vector<TopTree::triggeredObject> > trigfilters = event->getTriggerFilters();
            /////////////////////////////////
            // print when you change file  //
            /////////////////////////////////
            bool redotrigmap=false;
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename != currentFilename){
                previousFilename = currentFilename;
                iFile++;
                redotrigmap=true;
                cout<<"File changed!!! => iFile = "<<iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << dataSetName << endl;
            }
            
            // get run number
            int currentRun = event->runId();
            
            if(previousRun != currentRun){
                previousRun = currentRun;
                redotrigmap=true;
            }
            // get trigger info:
            
            
            for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermap.begin(); iter != triggermap.end(); iter++){
                if(redotrigmap){
		  Int_t loc= treeLoader.iTrigger(iter->first, currentRun,iFile);
                    iter->second.first=loc;
                }
                // and check if it fired:
                if(iter->second.first>=0 && iter->second.first!=9999) // trigger exists
                    iter->second.second=event->trigHLT(iter->second.first);
                else
                    iter->second.second=false;
            }
            
            // now check if the appropriate triggers fired for each analysis:
            trig_eplusjets=0;
            for(UInt_t itrig=0; itrig<ejetstriggers.size() && trig_eplusjets==0; itrig++){
                if(triggermap[ejetstriggers[itrig]].second)
                    trig_eplusjets=1;
            }
            trig_muplusjets=0;
            for(UInt_t itrig=0; itrig<mujetstriggers.size() && trig_muplusjets==0; itrig++){
                if(triggermap[mujetstriggers[itrig]].second)
                    trig_muplusjets=1;
            }
            trig_dilepton_ee=0;
            for(UInt_t itrig=0; itrig<eetriggers.size() && trig_dilepton_ee==0; itrig++){
                if(triggermap[eetriggers[itrig]].second)
                    trig_dilepton_ee=1;
            }
            trig_dilepton_emu=0;
            for(UInt_t itrig=0; itrig<emutriggers.size() && trig_dilepton_emu==0; itrig++){
                if(triggermap[emutriggers[itrig]].second)
                    trig_dilepton_emu=1;
            }
            trig_dilepton_mumu=0;
            for(UInt_t itrig=0; itrig<mumutriggers.size() && trig_dilepton_mumu==0; itrig++){
                if(triggermap[mumutriggers[itrig]].second)
                    trig_dilepton_mumu=1;
            }
            trig_displaced=0;
            for(UInt_t itrig=0; itrig<displtriggers.size() && trig_displaced==0; itrig++){
                if(triggermap[displtriggers[itrig]].second)
                    trig_displaced=1;
            }
            
            run_num=event->runId();
            evt_num=event->eventId();
            lumi_num=event->lumiBlockId();
            nvtx=vertex.size();
            npu=(int)event->nTruePU();
            
            if( run_num > 10000){//data
                isdata=1;
            }
            
            bookkeeping->Fill();
            
            
            
            if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) && !isdata ) {
                genjets = treeLoader.LoadGenJet(ievt,false);
                //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
            }
            else{
                isdata=1;
            }
            
            
            /////////////////////////////////
            // DETERMINE EVENT SCALEFACTOR //
            /////////////////////////////////
            
            // scale factor for the event
            float scaleFactor = 1.;
            
            // PU reweighting
            
            //            double lumiWeight = LumiWeights.ITweight( nvtx ); // simplest reweighting, just use reconstructed number of PV.
            double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() ); // reweighting using number of true pv
            
            if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
                lumiWeight=1;
            
            // filled into output file
            pu_weight=lumiWeight;
            
            scaleFactor = scaleFactor*lumiWeight;
            
            
            
            
            /////////////////////////////////////////////////////////////////////////////
            // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
            /////////////////////////////////////////////////////////////////////////////
            // not applied during CSA14/PHYS14
            //            if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
            //
            //                jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);
            
            /////////////////////
            // EVENT SELECTION //
            /////////////////////
            
            //Declare selection instance
            Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
            
            // the default selection is fine for normal use - if you want a special selection you can use the functions here
            //selection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); //  void setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon);
            //selection.setMuonCuts(20,2.5,1.0,2.0,0.3,1,0.5,5,0); // void setMuonCuts(float Pt, float Eta, float RelIso, float d0, float DRJets, int NMatchedStations, float Dz, int NTrackerLayersWithMeas, int NValidPixelHits);
            //            selection.setElectronCuts(20,2.5,1.0,2.0,0.5,0.4,0); // void setElectronCuts(float Pt, float Eta, float RelIso, float d0, float MVAId, float DRJets, int MaxMissingHits);
            
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
            
            if(!isGoodPV)
                continue;
            
            missingEt=mets[0]->Pt();
            
            // get the 'good' objects from the selection object
            vector<TRootPFJet*> selectedJets= selection.GetSelectedJets();
            vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons();
            vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons("Loose","PHYS14",true);
            vector<TRootElectron*> selectedMediumElectrons = selection.GetSelectedElectrons("Medium","PHYS14",true);
            vector<TRootElectron*> selectedTightElectrons = selection.GetSelectedElectrons("Tight","PHYS14",true);
            //            vector<TRootElectron*> selectedMediumElectronsMVA = selection.GetSelectedElectrons("Medium","PHYS14",false);
            //            vector<TRootElectron*> selectedTightElectronsMVA = selection.GetSelectedElectrons("Tight","PHYS14",false);
            
            
            vector<TRootElectron*> displacedElectrons = selection.GetSelectedDisplacedElectrons();
            vector<TRootMuon*> displacedMuons = selection.GetSelectedDisplacedMuons();
            // bookkeeping
            eleeffaverage[0]+=init_electrons.size()*scaleFactor;
            eleeffaverage[1]+=selectedElectrons.size()*scaleFactor;
            muoneffaverage[0]+=init_muons.size()*scaleFactor;
            muoneffaverage[1]+=selectedMuons.size()*scaleFactor;
            jeteffaverage[0]+=init_jets_corrected.size()*scaleFactor;
            jeteffaverage[1]+=selectedJets.size()*scaleFactor;
            
            
            // loop over electrons
            nElectrons=0;
            for(int iele=0; iele<selectedElectrons.size() && nElectrons<10; iele++){
                pT_electron[nElectrons]=selectedElectrons[iele]->Pt();
                phi_electron[nElectrons]=selectedElectrons[iele]->Phi();
                eta_electron[nElectrons]=selectedElectrons[iele]->Eta();
                E_electron[nElectrons]=selectedElectrons[iele]->E();
                d0_electron[nElectrons]=selectedElectrons[iele]->d0();
                loose_electron[nElectrons]=1;
                medium_electron[nElectrons]=0;
                tight_electron[nElectrons]=0;
                pfIso_electron[nElectrons]=selectedElectrons[iele]->relPfIso(3,0);
                charge_electron[nElectrons]=selectedElectrons[iele]->charge();
                for(int jele=0; jele<selectedMediumElectrons.size(); jele++){
                    if (selectedElectrons[iele]->DeltaR(*(selectedMediumElectrons[jele]))<0.001)
                        medium_electron[nElectrons]=1;
                }
                for(int jele=0; jele<selectedTightElectrons.size(); jele++){
                    if (selectedElectrons[iele]->DeltaR(*(selectedTightElectrons[jele]))<0.001)
                        tight_electron[nElectrons]=1;
                }
                //                for(int jele=0; jele<selectedMediumElectronsMVA.size(); jele++){
                //                    if (selectedElectrons[iele]->DeltaR(*(selectedMediumElectronsMVA[jele]))<0.001)
                //                        mediumMVA_electron[nElectrons]=1;
                //                }
                //                for(int jele=0; jele<selectedTightElectronsMVA.size(); jele++){
                //                    if (selectedElectrons[iele]->DeltaR(*(selectedTightElectronsMVA[jele]))<0.001)
                //                        tightMVA_electron[nElectrons]=1;
                //                }
                //
                
                nElectrons++;
            }
            // loop over muons
            nMuons=0;
            for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){
                pT_muon[nMuons]=selectedMuons[imuo]->Pt();
                phi_muon[nMuons]=selectedMuons[imuo]->Phi();
                eta_muon[nMuons]=selectedMuons[imuo]->Eta();
                E_muon[nMuons]=selectedMuons[imuo]->E();
                d0_muon[nMuons]=selectedMuons[imuo]->d0();
                pfIso_muon[nMuons]=selectedMuons[imuo]->relPfIso(3,0);
                
                
                charge_muon[nMuons]=selectedMuons[imuo]->charge();
                nMuons++;
            }
            // loop over jets
            nJets=0;
            for(int ijet=0; ijet<selectedJets.size() && nJets<20; ijet++){
                pT_jet[nJets]=selectedJets[ijet]->Pt();
                phi_jet[nJets]=selectedJets[ijet]->Phi();
                eta_jet[nJets]=selectedJets[ijet]->Eta();
                E_jet[nJets]=selectedJets[ijet]->E();
                dx_jet[nJets]=selectedJets[ijet]->vx();
                dy_jet[nJets]=selectedJets[ijet]->vy();
                btag_jet[nJets]=selectedJets[ijet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                nJets++;
            }
            
            
            if( nElectrons+nMuons>=2){
                myTree->Fill();
            }
            else if(nElectrons>=1&&nJets>3){
                myTree->Fill();
            }
            else if(nMuons >= 1 && nJets > 3){
                myTree->Fill();
            }
            // loop over electrons
            nElectrons=0;
            for(int iele=0; iele<displacedElectrons.size() && nElectrons<10; iele++){
                pT_electron[nElectrons]=displacedElectrons[iele]->Pt();
                phi_electron[nElectrons]=displacedElectrons[iele]->Phi();
                if(pT_electron[nElectrons]<20)
                     continue;
                eta_electron[nElectrons]=displacedElectrons[iele]->Eta();
                E_electron[nElectrons]=displacedElectrons[iele]->E();
                d0_electron[nElectrons]=displacedElectrons[iele]->d0();
                
                
                pfIso_electron[nElectrons]=displacedElectrons[iele]->relPfIso(3,0);
                charge_electron[nElectrons]=displacedElectrons[iele]->charge();
                nElectrons++;
            }
            // loop over muons
            nMuons=0;
            for(int imuo=0; imuo<displacedMuons.size() && nMuons<10; imuo++){
                pT_muon[nMuons]=displacedMuons[imuo]->Pt();
                phi_muon[nMuons]=displacedMuons[imuo]->Phi();
                if(pT_muon[nMuons]<20)
                    continue;
                eta_muon[nMuons]=displacedMuons[imuo]->Eta();
                E_muon[nMuons]=displacedMuons[imuo]->E();
                d0_muon[nMuons]=displacedMuons[imuo]->d0();
                pfIso_muon[nMuons]=displacedMuons[imuo]->relPfIso(3,0);
                
                
                charge_muon[nMuons]=displacedMuons[imuo]->charge();
                nMuons++;
            }
            // loop over jets
            nJets=0;
            for(int ijet=0; ijet<selectedJets.size() && nJets<20; ijet++){
                pT_jet[nJets]=selectedJets[ijet]->Pt();
                phi_jet[nJets]=selectedJets[ijet]->Phi();
                
                if(pT_jet[nJets]<20)
                    continue;
                
                eta_jet[nJets]=selectedJets[ijet]->Eta();
                E_jet[nJets]=selectedJets[ijet]->E();
                dx_jet[nJets]=selectedJets[ijet]->vx();
                dy_jet[nJets]=selectedJets[ijet]->vy();
                btag_jet[nJets]=selectedJets[ijet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                nJets++;
            }
            if(nElectrons+nMuons>=2){
                noselTree->Fill();
            }
            
        }			//loop on events
        
        cout << "number of (weighted) selected jets: " << jeteffaverage[1] ; if(jeteffaverage[0]>0) cout << " at efficiency of " << jeteffaverage[1]/jeteffaverage[0];
        cout << endl;
        cout << "number of (weighted) selected electrons: " << eleeffaverage[1] ; if(eleeffaverage[0]>0) cout << " at efficiency of " << eleeffaverage[1]/eleeffaverage[0];
        cout << endl;
        cout << "number of (weighted) selected muons: " << muoneffaverage[1] ; if(muoneffaverage[0]>0) cout << " at efficiency of " << muoneffaverage[1]/muoneffaverage[0];
        cout << endl;
        
        
        
        cout<<endl;
        
        
        
        //////////////
        // CLEANING //
        //////////////
        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        
        myTree->Write();
        fileout->Write();
        fileout->Close();
        //        delete myTree;
        delete fileout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }				//loop on datasets
    
    
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
