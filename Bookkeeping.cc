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
            
            run_num=event->runId();
            evt_num=event->eventId();
            lumi_num=event->lumiBlockId();
            nvtx=vertex.size();
            npu=(int)event->nTruePU();
            
            if( run_num > 10000){//data
                isdata=1;
            }
            
            bookkeeping->Fill();
            
	}        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        bookkeeping->Write();
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
    
    return 0;
}
