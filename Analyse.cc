///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Ntuple.cc: This macro is intended to be an example analysis macro which works out of the box.           /////
/////       It should serve as the first port of call for new users of the TopBrussels framework.             /////
/////      (in addition it is used by Freya for occasional studies when she has time)                         /////
/////     Last Modified: Mon 16 February 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <map>


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
	
	string fileList;
	string samplename="unknown";
    if (argc > 1)
        fileList = (string)argv[1];
	
	if (argc > 2)
		samplename= (string) argv[2];
	
	cout << "looping over files: " << fileList << ", and sample name: " << samplename << endl;
	float lumiweight=16.09+42.45;
	float crosssec=1;
	if(samplename.find("Data")!=string::npos){
		lumiweight=-1;
		
	}
	else{
		// MC samples... do other bookkeeping if needed.
		if(samplename.find("TTJets")!=string::npos){
			crosssec=831.76;
		}
		if(samplename.find("DY")!=string::npos){
				crosssec=2008.4*3;
		}
	}
		
	//////////////////////
	// straight from makeclass:
	/////////////////////
	
	Int_t           isdata;
	Int_t           run_num;
	Int_t           evt_num;
	Int_t           lumi_num;
	Int_t           nvtx;
	Int_t           npu;
	Int_t           trig_dilepton_emu;
	Int_t           trig_dilepton_ee;
	Int_t           trig_dilepton_mumu;
	Int_t           trig_eplusjets;
	Int_t           trig_muplusjets;
	Int_t           trig_displaced;
	Int_t           nElectrons;
	Double_t        pT_electron[5];   //[nElectrons]
	Double_t        phi_electron[5];   //[nElectrons]
	Double_t        eta_electron[5];   //[nElectrons]
	Double_t        E_electron[5];   //[nElectrons]
	Double_t        pfIso_electron[5];   //[nElectrons]
	Int_t           charge_electron[5];   //[nElectrons]
	Double_t        d0_electron[5];   //[nElectrons]
	Double_t        dz_electron[5];   //[nElectrons]
	Double_t        d0bs_electron[5];   //[nElectrons]
	Double_t        dzbs_electron[5];   //[nElectrons]
	Int_t           medium_electron[5];   //[nElectrons]
	Int_t           disp_electron[5];   //[nElectrons]
	Int_t           nMuons;
	Double_t        pT_muon[8];   //[nMuons]
	Double_t        phi_muon[8];   //[nMuons]
	Double_t        eta_muon[8];   //[nMuons]
	Double_t        E_muon[8];   //[nMuons]
	Double_t        pfIso_muon[8];   //[nMuons]
	Int_t           charge_muon[8];   //[nMuons]
	Int_t           disp_muon[8];   //[nMuons]
	Int_t           id_muon[8];   //[nMuons]
	Double_t        d0_muon[8];   //[nMuons]
	Double_t        dz_muon[8];   //[nMuons]
	Double_t        d0bs_muon[8];   //[nMuons]
	Double_t        dzbs_muon[8];   //[nMuons]
	Int_t           nJets;
	Double_t        pT_jet[6];   //[nJets]
	Double_t        phi_jet[6];   //[nJets]
	Double_t        eta_jet[6];   //[nJets]
	Double_t        E_jet[6];   //[nJets]
	Double_t        d0_jet[6];   //[nJets]
	Double_t        dz_jet[6];   //[nJets]
	Double_t        btag_jet[6];   //[nJets]
	Double_t        missingEt;
	Double_t        pu_weight;
	
	// List of branches
	TBranch        *b_isdata;   //!
	TBranch        *b_run_num;   //!
	TBranch        *b_evt_num;   //!
	TBranch        *b_lumi_num;   //!
	TBranch        *b_nvtx;   //!
	TBranch        *b_npu;   //!
	TBranch        *b_trig_dilepton_emu;   //!
	TBranch        *b_trig_dilepton_ee;   //!
	TBranch        *b_trig_dilepton_mumu;   //!
	TBranch        *b_trig_eplusjets;   //!
	TBranch        *b_trig_muplusjets;   //!
	TBranch        *b_trig_displaced;   //!
	TBranch        *b_nElectrons;   //!
	TBranch        *b_pT_electron;   //!
	TBranch        *b_phi_electron;   //!
	TBranch        *b_eta_electron;   //!
	TBranch        *b_E_electron;   //!
	TBranch        *b_pfIso_electron;   //!
	TBranch        *b_charge_electron;   //!
	TBranch        *b_d0_electron;   //!
	TBranch        *b_dz_electron;   //!
	TBranch        *b_d0bs_electron;   //!
	TBranch        *b_dzbs_electron;   //!
	TBranch        *b_medium_electron;   //!
	TBranch        *b_disp_electron;   //!
	TBranch        *b_nMuons;   //!
	TBranch        *b_pT_muon;   //!
	TBranch        *b_phi_muon;   //!
	TBranch        *b_eta_muon;   //!
	TBranch        *b_E_muon;   //!
	TBranch        *b_pfIso_muon;   //!
	TBranch        *b_charge_muon;   //!
	TBranch        *b_disp_muon;   //!
	TBranch        *b_id_muon;   //!
	TBranch        *b_d0_muon;   //!
	TBranch        *b_dz_muon;   //!
	TBranch        *b_d0bs_muon;   //!
	TBranch        *b_dzbs_muon;   //!
	TBranch        *b_nJets;   //!
	TBranch        *b_pT_jet;   //!
	TBranch        *b_phi_jet;   //!
	TBranch        *b_eta_jet;   //!
	TBranch        *b_E_jet;   //!
	TBranch        *b_d0_jet;   //!
	TBranch        *b_dz_jet;   //!
	TBranch        *b_btag_jet;   //!
	TBranch        *b_missingEt;   //!
	TBranch        *b_pu_weight;   //!
	
	
	TChain * ch = new TChain("tree","tree");
	ch->SetBranchAddress("isdata", &isdata, &b_isdata);
	ch->SetBranchAddress("run_num", &run_num, &b_run_num);
	ch->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
	ch->SetBranchAddress("lumi_num", &lumi_num, &b_lumi_num);
	ch->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
	ch->SetBranchAddress("npu", &npu, &b_npu);
	ch->SetBranchAddress("trig_dilepton_emu", &trig_dilepton_emu, &b_trig_dilepton_emu);
	ch->SetBranchAddress("trig_dilepton_ee", &trig_dilepton_ee, &b_trig_dilepton_ee);
	ch->SetBranchAddress("trig_dilepton_mumu", &trig_dilepton_mumu, &b_trig_dilepton_mumu);
	ch->SetBranchAddress("trig_eplusjets", &trig_eplusjets, &b_trig_eplusjets);
	ch->SetBranchAddress("trig_muplusjets", &trig_muplusjets, &b_trig_muplusjets);
	ch->SetBranchAddress("trig_displaced", &trig_displaced, &b_trig_displaced);
	ch->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
	ch->SetBranchAddress("pT_electron", pT_electron, &b_pT_electron);
	ch->SetBranchAddress("phi_electron", phi_electron, &b_phi_electron);
	ch->SetBranchAddress("eta_electron", eta_electron, &b_eta_electron);
	ch->SetBranchAddress("E_electron", E_electron, &b_E_electron);
	ch->SetBranchAddress("pfIso_electron", pfIso_electron, &b_pfIso_electron);
	ch->SetBranchAddress("charge_electron", charge_electron, &b_charge_electron);
	ch->SetBranchAddress("d0_electron", d0_electron, &b_d0_electron);
	ch->SetBranchAddress("dz_electron", dz_electron, &b_dz_electron);
	ch->SetBranchAddress("d0bs_electron", d0bs_electron, &b_d0bs_electron);
	ch->SetBranchAddress("dzbs_electron", dzbs_electron, &b_dzbs_electron);
	ch->SetBranchAddress("medium_electron", medium_electron, &b_medium_electron);
	ch->SetBranchAddress("disp_electron", disp_electron, &b_disp_electron);
	ch->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
	ch->SetBranchAddress("pT_muon", pT_muon, &b_pT_muon);
	ch->SetBranchAddress("phi_muon", phi_muon, &b_phi_muon);
	ch->SetBranchAddress("eta_muon", eta_muon, &b_eta_muon);
	ch->SetBranchAddress("E_muon", E_muon, &b_E_muon);
	ch->SetBranchAddress("pfIso_muon", pfIso_muon, &b_pfIso_muon);
	ch->SetBranchAddress("charge_muon", charge_muon, &b_charge_muon);
	ch->SetBranchAddress("disp_muon", disp_muon, &b_disp_muon);
	ch->SetBranchAddress("id_muon", id_muon, &b_id_muon);
	ch->SetBranchAddress("d0_muon", d0_muon, &b_d0_muon);
	ch->SetBranchAddress("dz_muon", dz_muon, &b_dz_muon);
	ch->SetBranchAddress("d0bs_muon", d0bs_muon, &b_d0bs_muon);
	ch->SetBranchAddress("dzbs_muon", dzbs_muon, &b_dzbs_muon);
	ch->SetBranchAddress("nJets", &nJets, &b_nJets);
	ch->SetBranchAddress("pT_jet", pT_jet, &b_pT_jet);
	ch->SetBranchAddress("phi_jet", phi_jet, &b_phi_jet);
	ch->SetBranchAddress("eta_jet", eta_jet, &b_eta_jet);
	ch->SetBranchAddress("E_jet", E_jet, &b_E_jet);
	ch->SetBranchAddress("d0_jet", d0_jet, &b_d0_jet);
	ch->SetBranchAddress("dz_jet", dz_jet, &b_dz_jet);
	ch->SetBranchAddress("btag_jet", btag_jet, &b_btag_jet);
	ch->SetBranchAddress("missingEt", &missingEt, &b_missingEt);
	ch->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
	/// end of stuff from MakeClass
	
	
	// now start doing stuff:
	ch->Add(fileList.c_str());
	
	TChain *bookkeeping = new TChain("startevents","startevents");
	bookkeeping->Add(fileList.c_str());
	
	Float_t normalisation=bookkeeping->GetEntries();
	
	Long64_t nentries = ch->GetEntries();
	std::cout << "looping over:" << nentries << " entries" << " with normalisation " << normalisation <<  ", cross section " << crosssec <<  endl;
	
	TFile *outfile = new TFile("histograms.root","update");
	outfile->mkdir(samplename.c_str());
	outfile->cd(samplename.c_str());
	// book histograms here!

	Int_t ncuts=5;
	TH2D *cutflow = new TH2D("cutflow","cutflow",ncuts,-0.5,ncuts-0.5,3,-0.5,2.5);
	cutflow->Sumw2();
	TH1D* electronIsoAll = new TH1D("electronIsoAll","pf Iso",100,0,1);
	electronIsoAll->Sumw2();
	TH1D* electronIsoPresel = (TH1D*) electronIsoAll->Clone("electronIsoPresel");
	TH1D* muonIsoAll = (TH1D*) electronIsoAll->Clone("muonIsoAll");
	TH1D* muonIsoPresel = (TH1D*) electronIsoAll->Clone("muonIsoPresel");

	TH1D* electrond0All = new TH1D("electrond0All","'|d0",100,-5,5);
	electrond0All->Sumw2();
	TH1D* electrond0Presel = (TH1D*) electrond0All->Clone("electrond0Presel");
	TH1D* muond0Presel = (TH1D*) electrond0All->Clone("muond0Presel");
	TH1D* muond0All = (TH1D*) electrond0All->Clone("muond0All");

	
	
	
	for(Long64_t jentry=0; jentry<normalisation; ++jentry){
		cutflow->Fill(0.,0.,1.);
		float evtweight=lumiweight;
		if(evtweight<0)
			evtweight=1;
		evtweight*=lumiweight;
		evtweight/=normalisation;
		cutflow->Fill(0.,1.,evtweight);
		cutflow->Fill(0.,2.,1./normalisation);
	}
	Long64_t nbytes = 0, nb = 0, ifile=0, ientry, iele, imuo, ijet, index_ele, index_mu;
	
	TLorentzVector worker1, worker2;
	float cutstep=1;
	for (Long64_t jentry=0; jentry<nentries;++jentry) {
		ientry = ch->LoadTree(jentry);
		if(ientry< 0)
			cout << "oops no tree loaded!" << endl;
		if(ifile!=ch->GetTreeNumber()){
			ifile=ch->GetTreeNumber();
		}
		ientry = ch->GetEntry(jentry);
		if(ientry< 0)
			cout << "oops no tree retrieved!" << endl;
		
		if(jentry%10000==0)
			cout << jentry << "/" << nentries << endl;
			
		double totalweight = 1;
		if(lumiweight>0){
			totalweight*=lumiweight;
			totalweight/=normalisation;
			totalweight*=pu_weight;
		}
		cutflow->Fill(cutstep,0.,1.);		cutflow->Fill(cutstep,1.,totalweight);		cutflow->Fill(cutstep,2.,1./normalisation); cutstep++;
		
//		cout << "number of :" << nElectrons << " " << jentry << "/" << nentries <<		" " << totalweight << endl;
		if(trig_displaced < 1 && lumiweight<0)// only do when data
			continue;
		
		cutflow->Fill(cutstep,0.,1.);		cutflow->Fill(cutstep,1.,totalweight);		cutflow->Fill(cutstep,2.,1./normalisation); cutstep++;
		
		if (nElectrons<1 || nMuons< 1)
			continue;
		
		cutflow->Fill(cutstep,0.,1.);		cutflow->Fill(cutstep,1.,totalweight);		cutflow->Fill(cutstep,2.,1./normalisation); cutstep++;
		
		index_ele=-1;
		for(iele=0; iele< nElectrons; iele++){
//			cout << pfIso_electron[iele] << "," << d0bs_electron[iele] << endl;

			electronIsoAll->Fill(pfIso_electron[iele],totalweight);
			electrond0All->Fill(d0bs_electron[iele],totalweight);
			
			if( disp_electron[iele]==0)
				continue;
			index_ele=iele;
			if(0){
			worker1.SetPtEtaPhiE(pT_electron[iele],eta_electron[iele],phi_electron[iele],E_electron[iele]);
			for(ijet=0; ijet<nJets; ijet++){
				// check for overlapping jets..
			}
			}
			
		}
		index_mu=-1;
		for(imuo=0; imuo< nMuons; imuo++){
//			cout << pfIso_muon[imuo] << "," << d0bs_muon[imuo] << endl;
			
			muonIsoAll->Fill(pfIso_muon[imuo],totalweight);
			muond0All->Fill(d0bs_muon[imuo],totalweight);
			
			if( disp_muon[imuo]==0)
				continue;
			
			index_mu=imuo;
			if(0){
				worker1.SetPtEtaPhiE(pT_muon[imuo],eta_muon[imuo],phi_muon[imuo],E_muon[imuo]);
				for(ijet=0; ijet<nJets; ijet++){
					// check for overlapping jets..
				}
			}
			
		}
		if( index_mu<0 || index_ele < 0)
			continue;
		
		cutflow->Fill(cutstep,0.,1.);		cutflow->Fill(cutstep,1.,totalweight);		cutflow->Fill(cutstep,2.,1./normalisation); cutstep++;
		
	}
	cout << "done, now closing" << endl;
	outfile->cd(samplename.c_str());
	TH1D *cutflow_noweights = (TH1D*) cutflow->ProfileY("cutflow_noweights",1,1);
	TH1D *cutflow_weights = (TH1D*) cutflow->ProfileY("cutflow_weights",2,2);
	TH1D *cutflow_eff = (TH1D*) cutflow->ProfileY("cutflow_eff",3,3);
	
	
	outfile->Write("kOverwrite");
	
    return 0;
}
