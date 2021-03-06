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
#include <map>

//used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "../TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"
#include "../TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TH2D.h"

using namespace std;
using namespace reweight;
using namespace TopTree;



Float_t calc_weight_p3(std::vector<float> effvals,std::vector<float> sfvals,int nobj, bool withsf){
    Float_t pdata= 0;
    if( effvals.size() < nobj ){
        return -1;
    }
    if( sfvals.size() < nobj ){
        return -1;
    }
    
    
    for(int ii=0; ii<nobj; ii++){
        Float_t workval1 = effvals[ii];
        if(withsf)
            workval1*=sfvals[ii];
        for(int jj=0; jj<nobj; jj++){
            if(ii==jj)
                continue;
            Float_t workval2 = effvals[jj];
            if(withsf)
                workval2*=sfvals[jj];
            for(int kk=0; kk<nobj; kk++){
                if (kk==ii)
                    continue;
                if (kk==jj)
                    continue;
                
                Float_t workval3 = effvals[kk];
                if(withsf)
                    workval3*=sfvals[kk];
                for(int ll=0; ll < nobj; ll++){
                    if (ll==kk)
                        continue;
                    if (ll==jj)
                        continue;
                    if (ll == ii)
                        continue;
                    
                    Float_t workval4 = effvals[ll];
                    if(withsf)
                        workval4*=sfvals[ll];
                    workval3*=(1-workval4);
                }
                pdata+=workval1*workval2*workval3;
            }
        }
    }
    return pdata;
}

Float_t calc_weight_p2(std::vector<float> effvals,std::vector<float> sfvals,int nobj, bool withsf){
    Float_t pdata= 0;
    if( effvals.size() < nobj ){
        return -1;
    }
    if( sfvals.size() < nobj ){
        return -1;
    }
    
    
    for(int ii=0; ii<nobj; ii++){
        Float_t workval1 = effvals[ii];
        if(withsf)
            workval1*=sfvals[ii];
        for(int jj=0; jj<nobj; jj++){
            if(ii==jj)
                continue;
            Float_t workval2 = effvals[jj];
            if(withsf)
                workval2*=sfvals[jj];
            for(int kk=0; kk<nobj; kk++){
                if (kk==ii)
                    continue;
                if (kk==jj)
                    continue;
                
                Float_t workval3 = effvals[kk];
                if(withsf)
                    workval3*=sfvals[kk];
                workval2*=(1.-workval3);
            }
            pdata+=workval1*workval2;
        }
        
    }
    return pdata;
}



Float_t calc_weight_p1(std::vector<float> effvals,std::vector<float> sfvals, int nobj, bool withsf){
    Float_t pdata= 0;
    if( effvals.size() < nobj ){
        return -1;
    }
    if( sfvals.size() < nobj ){
        return -1;
    }
    for(int ii=0; ii<nobj; ii++){
        Float_t workval1 = effvals[ii];
        if(withsf)
            workval1*=sfvals[ii];
        for(int jj=ii; jj<nobj; jj++){
            if(ii==jj)
                continue;
            Float_t workval2 = effvals[jj];
            if(withsf)
                workval2*=sfvals[jj];
            workval1*=1.-workval2;
        }
        pdata+=workval1;
        
    }
    return pdata;
}


Float_t calc_weight_p0(std::vector<float> effvals,std::vector<float> sfvals, int nobj, bool withsf){
    Float_t pdata= 1;
    if( effvals.size() < nobj ){
        return -1;
    }
    if( sfvals.size() < nobj ){
        return -1;
    }
    
    
    for(int ii=0; ii<nobj; ii++){
        //        cout << effvals[ii] << " " << sfvals[ii] << endl;
        Float_t workval=effvals[ii];
        if( withsf)
            workval *=sfvals[ii];
        pdata*= 1.-workval;
    }
    return pdata;
}




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
    cout << " - Load datasets ..." << endl;
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
    
    
    string pathToTTAB="/user/fblekman/DisplacedSusy2016/CMSSW_8_0_22/src/TopBrussels/HToZZBachelorProjectNtupleMaker/..";
    
    //Global variable
    //TRootEvent* event = 0;
    
    //nof selected events
    double NEvtsData = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
    ////////////////////////////////////
    //	Loop on datasets
    ////////////////////////////////////
    
    // last revisted on 11 dec 2015, https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
    std::vector<std::string> mujetstriggers;
    // newer triggers
    mujetstriggers.push_back("HLT_IsoMu18_v*");
    mujetstriggers.push_back("HLT_IsoMu17_eta2p1_v*");
    mujetstriggers.push_back("HLT_IsoMu20_v*");
    mujetstriggers.push_back("HLT_IsoTkMu20_v*");
    
    // older triggers
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_v*");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v*");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_v*");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*");
    mujetstriggers.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v*");
    
    std::vector<std::string> ejetstriggers;
    ejetstriggers.push_back("HLT_Ele23_WPLoose_Gsf_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v*");
    ejetstriggers.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v*");
    ejetstriggers.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v*");
    
    std::vector<std::string> eetriggers;
    eetriggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    eetriggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    
    std::vector<std::string> mumutriggers;
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
    mumutriggers.push_back("HLT_IsoMu20_eta2p1_v*");
    mumutriggers.push_back("HLT_IsoMu20_eta2p1_v*");
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    mumutriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
    
    std::vector<std::string> emutriggers;
    emutriggers.push_back("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
    emutriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*");
    emutriggers.push_back("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
    emutriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*");
    
    std::vector<std::string> displtriggers;
    displtriggers.push_back("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL_v*");
    displtriggers.push_back("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL_v*");
    
    std::vector<std::string> alltriggers;
    
    std::map<std::string,std::pair<int,bool> > triggermap;

    for(UInt_t itrig=0; itrig<mujetstriggers.size(); itrig++){
        triggermap[mujetstriggers[itrig]]=std::pair<int,bool>(-999,false);
        alltriggers.push_back(mujetstriggers[itrig]);
        
    }
    
    
    // all sorts of calibration loading:
    BTagCalibration btagsfweight_hf("CSVv2",pathToTTAB+"/TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_80X_ichep_incl_ChangedTo_mujets.csv");
    BTagCalibrationReader btagsfreader(&btagsfweight_hf,BTagEntry::OP_MEDIUM,"mujets","central");
    BTagCalibrationReader btagsfreadercomb(&btagsfweight_hf,BTagEntry::OP_MEDIUM,"mujets","central");
    
    
    MuonSFWeight musfweight(pathToTTAB+"/TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonID_Z_RunBCD_prompt80X_7p65.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");
    MuonSFWeight musfweightIso(pathToTTAB+"/TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonIso_Z_RunBCD_prompt80X_7p65.root","MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");

    ElectronSFWeight elesfweight(pathToTTAB+"/TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/egammaEffi.txt_SF2D_CutBasedMediumID.root","EGamma_SF2D");
    
    cout << "initiating lumi weights: " << endl;
    LumiReWeighting LumiWeights(pathToTTAB+"/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring16MiniAODv2-Asympt.root",pathToTTAB+"/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root","pileup","pileup");
    cout << "done loading lumi weights " << endl;
    //  TFile *file_pu_data = new TFile("/user/fblekman/TopTreeWork76X/TTP76Xv1/CMSSW_7_6_5/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert_Silver.root","READ");
    //  file_pu_data->cd();
    //  TH1D* pu_weight_histo = (TH1D*) file_pu_data->Get("pileup");
    //  pu_weight_histo->SetDirectory(0);
    //  pu_weight_histo->Scale(1. / pu_weight_histo->Integral());
    //  pu_weight_histo->SetName("pu_weight_histo");
    //  TFile *file_pu_mc = new TFile("/user/fblekman/TopTreeWork76X/TTP76Xv1/CMSSW_7_6_5/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root","READ");
    //  file_pu_mc->cd();
    //  TH1D* pu_weight_histo_denom = (TH1D*) file_pu_mc->Get("pileup");
    //  pu_weight_histo_denom->Scale(1./ pu_weight_histo_denom->Integral());
    //  pu_weight_histo->Divide(pu_weight_histo_denom);
    //  cout << "PU INFO, MANUALLY CALCULATED: " << endl;
    //  for(int ibin=1; ibin<pu_weight_histo->GetNbinsX(); ibin++)
    //    cout <<"bin " << ibin << " " << pu_weight_histo->GetBinContent(ibin) << endl;
    
    
    
    TFile *btagcalibs = new TFile("/user/fblekman/TopTreeWork76X/TTP76Xv1/CMSSW_7_6_5/src/TopBrussels/HToZZBachelorProjectNtupleMaker/output_btagging.root","READ");
    btagcalibs->cd();
    
    std::map<int,TH2D*> btag_efficiencies;
    btag_efficiencies[0]=(TH2D*)btagcalibs->Get("b_jetPtEta_Medium_LF");
    btag_efficiencies[4]=(TH2D*)btagcalibs->Get("b_jetPtEta_Medium_C");
    btag_efficiencies[5]=(TH2D*)btagcalibs->Get("b_jetPtEta_Medium_HF");
    btag_efficiencies[10]=(TH2D*)btagcalibs->Get("b_jetPtEta_LF");
    btag_efficiencies[14]=(TH2D*)btagcalibs->Get("b_jetPtEta_C");
    btag_efficiencies[15]=(TH2D*)btagcalibs->Get("b_jetPtEta_HF");
    
    
    
    
    
    std::vector<float> btag_vals(20,0);
    std::vector<float> btag_sf_vals(20,0);
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
        Int_t id_electron[10];
        Double_t d0bs_electron[10];
        Double_t dzbs_electron[10];
        Double_t pfIso_electron[10];
        Double_t sf_electron[10];
        Double_t sfh_electron[10];
        Double_t sfl_electron[10];
        Double_t drJet_electron[10];
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
        Double_t d0bs_muon[10];
        Double_t dzbs_muon[10];
        Double_t pfIso_muon[10];
        Int_t id_muon[10];
        Int_t loose_muon[10];
        Int_t medium_muon[10];
        Int_t tight_muon[10];
        
        Double_t sf_muon[10];
        Double_t sfh_muon[10];
        Double_t sfl_muon[10];
        Double_t sfIso_muon[10];
        Double_t sfIsoh_muon[10];
        Double_t sfIsol_muon[10];
        Double_t drJet_muon[10];
        Int_t charge_muon[10];
        
        Int_t nJets;
        Double_t pT_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Double_t dx_jet[20];
        Double_t dy_jet[20];
        Double_t sfb_jet[20];
        Double_t effb_jet[20];
        Double_t btag_jet[20];
        Int_t flav_jet[20];
        Double_t missingEt;
        // various weights
        Double_t pu_weight;
        Double_t btv_weight;
        Double_t lep_weight;
        Double_t trig_weight;
        Double_t mc_baseweight;
        Double_t mc_scaleupweight;
        Double_t mc_scaledownweight;
        Double_t mc_muonsfweight[3];
        Double_t mc_muonsfweightIso[3];
        Double_t mc_elesfweight[3];
        Double_t mc_btgsfweight1[5];
        Double_t mc_btgsfweight2[5];
        Double_t mc_btgsfweight3[5];
        Double_t mc_btgsfweight4[5];
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
        Int_t triggers_container[200];
        
        // cutflow variables:
        Int_t cutflow_pvcut;
        Int_t cutflow_ngoodele;
        Int_t cutflow_ngoodmuo;
        Int_t cutflow_ngoodjet;
        Int_t cutflow_ngoodeledisp;
        Int_t cutflow_ngoodjetdisp;
        Int_t cutflow_ngoodmuodisp;
        
        TTree *bookkeeping = new TTree("startevents","startevents");
        bookkeeping->Branch("run_num",&run_num,"run_num/I");
        bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
        bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
        bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
        bookkeeping->Branch("npu",&npu,"npu/I");
        bookkeeping->Branch("mc_baseweight",&mc_baseweight,"mc_baseweight/D");
        bookkeeping->Branch("cutflow_pvcut",&cutflow_pvcut,"cutflow_pvcut/I");
        bookkeeping->Branch("cutflow_ngoodele",&cutflow_ngoodele,"cutflow_ngoodele/I");
        bookkeeping->Branch("cutflow_ngoodmuo",&cutflow_ngoodmuo,"cutflow_ngoodmuo/I");
        bookkeeping->Branch("cutflow_ngoodmuodisp",&cutflow_ngoodmuodisp,"cutflow_ngoodmuodisp/I");
        bookkeeping->Branch("cutflow_ngoodeledisp",&cutflow_ngoodeledisp,"cutflow_ngoodeledisp/I");
        bookkeeping->Branch("cutflow_ngoodjet",&cutflow_ngoodjet,"cutflow_ngoodjet/I");
        bookkeeping->Branch("cutflow_ngoodjetdisp",&cutflow_ngoodjetdisp,"cutflow_ngoodjetdisp/I");
        
        
        
        
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
        myTree->Branch("drJet_electron",drJet_electron,"drJet_electron[nElectrons]/D");
        myTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
        myTree->Branch("sfl_electron",sfl_electron,"sfl_electron[nElectrons]/D");
        myTree->Branch("sfh_electron",sfh_electron,"sfh_electron[nElectrons]/D");
        myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        myTree->Branch("d0bs_electron",d0bs_electron,"d0b_electron[nElectrons]/D");
        myTree->Branch("dzbs_electron",dzbs_electron,"dzbs_electron[nElectrons]/D");
        myTree->Branch("loose_electron",loose_electron,"loose_electron[nElectrons]/I");
        myTree->Branch("medium_electron",medium_electron,"medium_electron[nElectrons]/I");
        myTree->Branch("tight_electron",tight_electron,"tight_electron[nElectrons]/I");
        //        myTree->Branch("mediumMVA_electron",mediumMVA_electron,"mediumMVA_electron[nElectrons]/I");
        //        myTree->Branch("tightMVA_electron",tightMVA_electron,"tightMVA_electron[nElectrons]/I");
        //
        
        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        myTree->Branch("sfl_muon",sfl_muon,"sfl_muon[nMuons]/D");
        myTree->Branch("sfh_muon",sfh_muon,"sfh_muon[nMuons]/D");
        myTree->Branch("sfIso_muon",sfIso_muon,"sfIso_muon[nMuons]/D");
        myTree->Branch("sfIsol_muon",sfIsol_muon,"sfIsol_muon[nMuons]/D");
        myTree->Branch("sfIsoh_muon",sfIsoh_muon,"sfIsoh_muon[nMuons]/D");
        myTree->Branch("pT_muon",pT_muon,"pT_muon[nMuons]/D");
        myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        myTree->Branch("drJet_muon",drJet_muon,"drJet_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        myTree->Branch("dzbs_muon",dzbs_muon,"dzbs_muon[nMuons]/D");
        myTree->Branch("d0bs_muon",d0bs_muon,"d0bs_muon[nMuons]/D");
        
        
        myTree->Branch("nJets",&nJets, "nJets/I");
        myTree->Branch("pT_jet",pT_jet,"pT_jet[nJets]/D");
        myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        myTree->Branch("sfb_jet",sfb_jet,"sfb_jet[nJets]/D");
        myTree->Branch("effb_jet",effb_jet,"effb_jet[nJets]/D");
        myTree->Branch("dx_jet",dx_jet,"dx_jet[nJets]/D");
        myTree->Branch("dy_jet",dy_jet,"dy_jet[nJets]/D");
        myTree->Branch("btag_jet",btag_jet,"btag_jet[nJets]/D");
        myTree->Branch("flav_jet",flav_jet,"flav_jet[nJets]/I");
        
        myTree->Branch("missingEt",&missingEt,"missingEt/D");
        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
        myTree->Branch("mc_baseweight",&mc_baseweight,"mc_baseweight/D");
        myTree->Branch("mc_scaleupweight",&mc_scaleupweight,"mc_scaleupweight/D");
        myTree->Branch("mc_scaledownweight",&mc_scaledownweight,"mc_scaledownweight/D");
        //        myTree->Branch("mc_muonsfweight",mc_muonsfweight,"mc_muonsfweight[3]/D");
        //        myTree->Branch("mc_elesfweight",mc_elesfweight,"mc_elesfweight[3]/D");
        myTree->Branch("mc_btgsfweight1",mc_btgsfweight1,"mc_btgsfweight1[5]/D");
        myTree->Branch("mc_btgsfweight2",mc_btgsfweight2,"mc_btgsfweight2[5]/D");
        myTree->Branch("mc_btgsfweight3",mc_btgsfweight3,"mc_btgsfweight3[5]/D");
        myTree->Branch("mc_btgsfweight4",mc_btgsfweight4,"mc_btgsfweight4[5]/D");
        
        
        myTree->Branch("id_electron",id_electron,"id_electron[nElectrons]/I");
        myTree->Branch("id_muon",id_muon,"id_muon[nMuons]/I");
        
                
        // automated branch adding for individual triggers:
        // limited by size of triggers_container!!
        cout << "checking triggers: " << endl;
        for(int iter_trig=0; iter_trig<alltriggers.size() && iter_trig<200; iter_trig++){
            TString trigname = alltriggers[iter_trig];
            trigname.ReplaceAll("_v*","");
            TString branchname = trigname+"/I";
            std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
            myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
            
        }
        
        cout << "done checking triggers: " << endl;
        
        
        
        Double_t workleptoneta, workleptonpt;
        
        //open files and load
        treeLoader.LoadDataset (datasets[d], anaEnv);
        
        
        //////////////////////////////////////////////////
        /// Initialize Jet energy correction factors   ///
        //////////////////////////////////////////////////
        
        vector<JetCorrectorParameters> vCorrParam;
        
        
        // data: /user/fblekman/DisplacedSusy2016/CMSSW_8_0_22/src/TopBrussels/HToZZBachelorProjectNtupleMaker/../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_UncertaintySources_AK4PFchs.txt
        // MC: /user/fblekman/DisplacedSusy2016/CMSSW_8_0_22/src/TopBrussels/HToZZBachelorProjectNtupleMaker/../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_Uncertainty_AK4PFchs.txt
//        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(pathToTTAB+"/TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_Uncertainty_AK4PFchs.txt", "Total")));
//        
//        // true means redo also the L1 corrections (see CMS documentation to learn what this means)
//        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
//        
//        
        ////////////////////////////////////
        //	Loop on events
        ////////////////////////////////////
        
        // some bookkeeping variables
        nEvents[d] = 0;
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
        
        // some printout
        cout << "running over " << datasets[d]->NofEvtsToRunOver() << endl;
        
        // run information
        TRootRun *runInfos = new TRootRun();
        datasets[d]->runTree()->SetBranchStatus("runInfos*",1);
        datasets[d]->runTree()->SetBranchAddress("runInfos",&runInfos);
        
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
            
            if(ievt%1000 == 0 || ievt < 10)
                std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected" << flush<<"\r";
            
            ////////////////
            // LOAD EVENT //
            ////////////////
            
	    //            cout << "loading event " << ievt << " !" << endl;
            TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
	    //            cout << "done loading event " << ievt << "! " << endl;
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
                datasets[d]->runTree()->GetEntry(iFile);
                
            }
            
            // get run number
            int currentRun = event->runId();
            
            if(previousRun != currentRun){
                previousRun = currentRun;
                redotrigmap=true;
                datasets[d]->runTree()->GetEntry(iFile);
            }
            
            
            if(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 1")>=0){
                mc_baseweight = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 1")))/(abs(event->originalXWGTUP()));
                mc_scaleupweight = event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 5"))/(abs(event->originalXWGTUP()));
                mc_scaledownweight = event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 9"))/(abs(event->originalXWGTUP()));
            }
            else if(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 1")>=0){
                mc_baseweight = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 1")))/(abs(event->originalXWGTUP()));
                mc_scaleupweight = event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 5"))/(abs(event->originalXWGTUP()));
                mc_scaledownweight = event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 9"))/(abs(event->originalXWGTUP()));
                
            }
            else{
                mc_baseweight=mc_scaleupweight=mc_scaledownweight=1.;
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
            
            for(UInt_t itrig=0; itrig<alltriggers.size() && itrig<200; itrig++){
                triggers_container[itrig]=0;
                if(triggermap[alltriggers[itrig]].second)
                    triggers_container[itrig]=1;
            }
            
            run_num=event->runId();
            evt_num=event->eventId();
            lumi_num=event->lumiBlockId();
            nvtx=vertex.size();
            npu=(int)event->nTruePU();
            
            if( run_num > 10000){//data
                isdata=1;
            }
            
            
            
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
            double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() ); // reweighting using number of true pv - unexpected behaviour so now just using histograms.
            //	double lumiWeight = pu_weight_histo->GetBinContent(pu_weight_histo->FindBin( (int)event->nTruePU() ));
            
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
            
            
            
            missingEt=mets[0]->Pt();
            
            // get the 'good' objects from the selection object
            vector<TRootPFJet*> selectedJets= selection.GetSelectedJets();
            vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons();
            vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons("Loose","Spring16_80X",true);
            vector<TRootElectron*> selectedMediumElectrons = selection.GetSelectedElectrons("Medium","Spring16_80X",true);
            vector<TRootElectron*> selectedTightElectrons = selection.GetSelectedElectrons("Tight","Spring16_80X",true);
            //            vector<TRootElectron*> selectedMediumElectronsMVA = selection.GetSelectedElectrons("Medium","PHYS14",false);
            //            vector<TRootElectron*> selectedTightElectronsMVA = selection.GetSelectedElectrons("Tight","PHYS14",false);
            
            
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
            cutflow_pvcut=isGoodPV;
            
            
            // loop over electrons
            nElectrons=0;
            for(int iele=0; iele<selectedElectrons.size() && nElectrons<10; iele++){
                if ( selectedElectrons[iele]->Pt() < 20)
                    continue;
                drJet_electron[nElectrons]=20;
                pT_electron[nElectrons]=selectedElectrons[iele]->Pt();
                phi_electron[nElectrons]=selectedElectrons[iele]->Phi();
                eta_electron[nElectrons]=selectedElectrons[iele]->Eta();
                E_electron[nElectrons]=selectedElectrons[iele]->E();
                d0_electron[nElectrons]=selectedElectrons[iele]->d0();
                d0bs_electron[nElectrons]=selectedElectrons[iele]->d0BeamSpot();
                //                dzbs_electron[nElectrons]=selectedElectrons[iele]->dzBeamSpot();
                loose_electron[nElectrons]=1;
                medium_electron[nElectrons]=0;
                tight_electron[nElectrons]=0;
                workleptoneta=selectedElectrons[iele]->Eta();
                workleptonpt=selectedElectrons[iele]->Pt();
                if(selectedElectrons[iele]->Pt()>500)
                    workleptonpt=499.999;
                else if(selectedElectrons[iele]->Pt()<20)
                    workleptonpt=20.0001;
                if(fabs(selectedElectrons[iele]->Eta())>2.4)
                    workleptoneta=2.3999;
                sf_electron[nElectrons]=elesfweight.at(workleptoneta,workleptonpt);
                sfh_electron[nElectrons]=elesfweight.at(workleptoneta,workleptonpt,1);
                sfl_electron[nElectrons]=elesfweight.at(workleptoneta,workleptonpt,-1);
                
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
            // get SFs for electrons
            // loop over muons
            nMuons=0;
            for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){
                if (selectedMuons[imuo]->Pt() < 20)
                    continue;
                drJet_muon[nMuons]=20;
                pT_muon[nMuons]=selectedMuons[imuo]->Pt();
                phi_muon[nMuons]=selectedMuons[imuo]->Phi();
                eta_muon[nMuons]=selectedMuons[imuo]->Eta();
                E_muon[nMuons]=selectedMuons[imuo]->E();
                d0_muon[nMuons]=selectedMuons[imuo]->d0();
                d0bs_muon[nMuons]=selectedMuons[imuo]->d0BeamSpot();
                //                dzbs_muon[nMuons]=selectedMuons[imuo]->dzBeamSpot();
                //                pfIso_muon[nMuons]=selectedMuons[imuo]->relPfIso(4,0);
                pfIso_muon[nMuons]=(selectedMuons[imuo]->chargedHadronIso(4) + max( 0.0, selectedMuons[imuo]->neutralHadronIso(4) + selectedMuons[imuo]->photonIso(4) - 0.5*selectedMuons[imuo]->puChargedHadronIso(4) ) ) / selectedMuons[imuo]->Pt(); // dBeta corrected
                workleptonpt=selectedMuons[imuo]->Pt();
                workleptoneta=selectedMuons[imuo]->Eta();
                if(selectedMuons[imuo]->Pt()>120)
                    workleptonpt=119.999;
                else if(selectedMuons[imuo]->Pt()<20)
                    workleptonpt=20.0001;
                
                if(workleptoneta>2.4)
                    workleptoneta=2.39999;
                if(workleptoneta<-2.4)
                    workleptoneta=-2.39999;
                
                sf_muon[nMuons]=musfweight.at(workleptoneta,workleptonpt);
                sfh_muon[nMuons]=musfweight.at(workleptoneta,workleptonpt,1);
                sfl_muon[nMuons]=musfweight.at(workleptoneta,workleptonpt,-1);
                sfIso_muon[nMuons]=musfweightIso.at(workleptoneta,workleptonpt);
                sfIsoh_muon[nMuons]=musfweightIso.at(workleptoneta,workleptonpt,1);
                sfIsol_muon[nMuons]=musfweightIso.at(workleptoneta,workleptonpt,-1);
                
                
                
                charge_muon[nMuons]=selectedMuons[imuo]->charge();
                nMuons++;
            }
            // loop over jets
            nJets=0;
            for(int ijet=0; ijet<selectedJets.size() && nJets<20; ijet++){
                bool badjet=false;
                for( int imuo=0; imuo<nMuons; imuo++){
                    worker.SetPtEtaPhiE(pT_muon[imuo],eta_muon[imuo],phi_muon[imuo],E_muon[imuo]);
                    if(selectedJets[ijet]->DeltaR(worker)<drJet_muon[imuo])
                        drJet_muon[imuo]=selectedJets[ijet]->DeltaR(worker);
                    if(selectedJets[ijet]->DeltaR(worker)<0.4)
                        badjet=true;
                }
                for( int iele=0; iele<nElectrons; iele++){
                    worker.SetPtEtaPhiE(pT_electron[iele],eta_electron[iele],phi_electron[iele],E_electron[iele]);
                    if(selectedJets[ijet]->DeltaR(worker)<drJet_electron[iele])
                        drJet_electron[iele]=selectedJets[ijet]->DeltaR(worker);
                    if(selectedJets[ijet]->DeltaR(worker)<0.4)
                        badjet=true;
                }
                if(badjet)
                    continue;
                pT_jet[nJets]=selectedJets[ijet]->Pt();
                phi_jet[nJets]=selectedJets[ijet]->Phi();
                eta_jet[nJets]=selectedJets[ijet]->Eta();
                E_jet[nJets]=selectedJets[ijet]->E();
                dx_jet[nJets]=selectedJets[ijet]->vx();
                dy_jet[nJets]=selectedJets[ijet]->vy();
                btag_jet[nJets]=selectedJets[ijet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                flav_jet[nJets]=selectedJets[ijet]->hadronFlavour();
                sfb_jet[nJets]=1.0;
                //cout << "jet flavour and SF: " << flav_jet[nJets] << endl;
                if(fabs(eta_jet[nJets])>2.4)
                    sfb_jet[nJets]=0.0;
                else{
                    if(fabs(flav_jet[nJets])==5)
                        sfb_jet[nJets]=btagsfreader.eval(selectedJets[ijet]);
                    else if(fabs(flav_jet[nJets])==4)
                        sfb_jet[nJets]=btagsfreader.eval(selectedJets[ijet]);
                    else if(fabs(flav_jet[nJets])==0)
                        sfb_jet[nJets]=btagsfreadercomb.eval(selectedJets[ijet]);
                }
                
                int binx =btag_efficiencies[flav_jet[nJets]]->GetXaxis()->FindBin(pT_jet[nJets]);
                int biny =btag_efficiencies[flav_jet[nJets]]->GetYaxis()->FindBin(eta_jet[nJets]);
                
                effb_jet[nJets]=btag_efficiencies[flav_jet[nJets]]->GetBinContent(binx,biny);
                if(btag_efficiencies[flav_jet[nJets]+10]->GetBinContent(binx,biny)>0)
                    effb_jet[nJets]/=btag_efficiencies[flav_jet[nJets]+10]->GetBinContent(binx,biny);
                
                
                
                nJets++;
                
            }
            // copy over to vectors and calculate SFs for btagging, inefficient but needs copy anyway so only using good ejts
            for(int ijet=0; ijet<nJets; ijet++){
                btag_sf_vals[ijet]=sfb_jet[ijet];
                btag_vals[ijet]=effb_jet[ijet];
                
            }
            float p0=calc_weight_p0(btag_vals,btag_sf_vals,nJets,false);
            float p0_sf = calc_weight_p0(btag_vals,btag_sf_vals,nJets,true);
            float p1=calc_weight_p1(btag_vals,btag_sf_vals,nJets,false);
            float p1_sf = calc_weight_p1(btag_vals,btag_sf_vals,nJets,true);
            float p2=calc_weight_p2(btag_vals,btag_sf_vals,nJets,false);
            float p2_sf = calc_weight_p2(btag_vals,btag_sf_vals,nJets,true);
            float p3=calc_weight_p3(btag_vals,btag_sf_vals,nJets,false);
            float p3_sf = calc_weight_p3(btag_vals,btag_sf_vals,nJets,true);
            if(1.-p0>0.)
                mc_btgsfweight1[0]=(1.-p0_sf)/(1.-p0);
            else
                mc_btgsfweight1[0]=0;
            if(1.-p0-p1>0.)
                mc_btgsfweight2[0]=(1.-p0_sf-p1_sf)/(1.-p0-p1);
            else
                mc_btgsfweight2[0]=0;
            
            if(1.-p0-p1-p2>0.)
                mc_btgsfweight3[0]=(1.-p0_sf-p1_sf-p2_sf)/(1.-p0-p1-p2);
            else
                mc_btgsfweight3[0]=0;
            
            if(1.-p0-p1-p2-p3>0.)
                mc_btgsfweight4[0]=(1.-p0_sf-p1_sf-p2_sf-p3_sf)/(1.-p0-p1-p2-p3);
            else
                mc_btgsfweight4[0]=0;
            
            //            cout << "values btag p0: " << p0 << ", " << p0_sf <<" " << p1 << " " << p1_sf   << " " << p2 << " " << p2_sf   <<" " << mc_btgsfweight1[0] << " " << mc_btgsfweight2[0] <<  endl;
            
            
            if( isGoodPV){
                
                if(nJets>=7 && nMuons<2 && nElectrons<2){
                    myTree->Fill();
                }
            }
            bookkeeping->Fill();
            
	    //            cout << "done with event " << ievt << endl;
            
        }//loop on events
        
        
        
        //////////////
        // CLEANING //
        //////////////
//        
//        if (jecUnc) delete jecUnc;
//        if (jetTools) delete jetTools;
        myTree->Write();
        fileout->Write();
        fileout->Close();
        //        delete myTree;
        delete fileout;
        
        //important: free memory    treeLoader.UnLoadDataset();
        
    }				//loop on datasets
    
    
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    btagcalibs->Close();
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
