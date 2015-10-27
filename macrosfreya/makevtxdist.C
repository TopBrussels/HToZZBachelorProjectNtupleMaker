{
	TH1D *pileup = new TH1D("pileup","pileup",60,0,60);
	pileup->Sumw2();
	TH1D *nvtxdist = (TH1D*) pileup->Clone("nvtxdist");
	
	TChain *ch = new TChain("startevents","startevents");

       
	// ch->Add("../DoubleEG-Run2015D-PromptReco-v3-CMSSW_74X_v8-74X_dataRun2_Prompt_v2_Ntupler/*.root");
	// ch->Add("../DoubleEG-Run2015D-PromptReco-v4-CMSSW_74X_v8-74X_dataRun2_Prompt_v2_Ntupler/*.root");
	// ch->Add("../DoubleMuon-Run2015D-PromptReco-v3-CMSSW_74X_v8-74X_dataRun2_Prompt_v2_Ntupler/*.root");
	// ch->Add("../DoubleMuon-Run2015D-PromptReco-v4-CMSSW_74X_v8-74X_dataRun2_Prompt_v2_Ntupler/*.root");
	// ch->Add("../MuonEG-Run2015D-PromptReco-v3-CMSSW_74X_v8-74X_dataRun2_Prompt_v2_Ntupler/*.root");
	// ch->Add("../MuonEG-Run2015D-PromptReco-v4-CMSSW_74X_v8-74X_dataRun2_Prompt_v2_Ntupler/*.root");
	// ch->Add("../SingleElectron-Run2015D-PromptReco-v3_Ntupler/*.root");
	// ch->Add("../SingleElectron-Run2015D-PromptReco-v4_Ntupler/*.root");
	// ch->Add("../SingleMuon-Run2015D-PromptReco-v3_Ntupler/*.root");
	// ch->Add("../SingleMuon-Run2015D-PromptReco-v4_Ntupler/*.root");


	
       	ch->Add("../DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Ntupler/*.root");
	ch->Add("../TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Ntupler/*.root");
	ch->Add("../WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Ntupler/*.root");
	
	std::cout << "using " << ch->GetEntries() << " events to make pu histogram:" << std::endl;
	ch->Draw("nvtx>>pileup"); 
	ch->Draw("nvtx>>nvtxdist");
//	ch->Print();
	
	double scaler = 1.;
	if(ch->GetEntries()>0)
		scaler /= ch->GetEntries();
		
	pileup->Scale(scaler);
	pileup->SaveAs("pileup_MC_RunIISpring15DR74-Asympt25ns.root");
	//       	pileup->SaveAs("nominal.root");
		
}
