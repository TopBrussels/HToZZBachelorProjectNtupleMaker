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
	
	
	for (unsigned int d = 0; d < datasets.size (); d++) {
		
		string previousFilename = "";
		int iFile = -1;
		string dataSetName = datasets[d]->Name();
		
		
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
		
		
		// various weights
		Double_t pu_weight;
		Int_t run_num;
		Int_t evt_num;
		Int_t lumi_num;
		Int_t nvtx;
		Int_t npu;
		Int_t isdata;
		
		TTree *bookkeeping = new TTree("startevents","startevents");
		bookkeeping->Branch("run_num",&run_num,"run_num/I");
		bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
		bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
		bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
		bookkeeping->Branch("npu",&npu,"npu/I");
		bookkeeping->Branch("isdata",&isdata,"isdata/I");
		
		
		//open files and load
		treeLoader.LoadDataset (datasets[d], anaEnv);
		
		
		//////////////////////////////////////////////////
		/// Initialize Jet energy correction factors   ///
		//////////////////////////////////////////////////
		
		////////////////////////////////////
		//	Loop on events
		////////////////////////////////////
		
		// some bookkeeping variables
		nEvents[d] = 0;
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
			
			bool redotrigmap=false;
			string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
		
			
			// get run number
			run_num=event->runId();
			evt_num=event->eventId();
			lumi_num=event->lumiBlockId();
			nvtx=vertex.size();
			npu=(int)event->nTruePU();
			isdata=0;
			if( run_num > 10000){//data
				isdata=1;
			}
			
			bookkeeping->Fill();
			
		}
		
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
