{
	TH1F *pileup = new TH1F("pileup","pileup",100,0,100);
	pileup->Sumw2();
	TH1F *nvtxdist = (TH1F*) pileup->Clone("nvtxdist");
	
	TChain *ch = new TChain("startevents","startevents");
//	ch->Add("../SingleMuon-Run2015C-23Sep2015-v1/*.root");
	ch->Add("../SingleMuon-Run2015D-PromptReco-v3/*.root");
	ch->Add("../SingleMuon-Run2015D-PromptReco-v4/*.root");
	
	ch->Add("../DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*.root");
	ch->Add("../TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*.root");
	
	std::cout << "using " << ch->GetEntries() << " events to make pu histogram:" << std::endl;
	ch->Draw("npu>>pileup");
	ch->Draw("nvtx>>nvtxdist");
//	ch->Print();
	
	double scaler = 1.;
	if(ch->GetEntries()>0)
		scaler /= ch->GetEntries();
		
	pileup->Scale(scaler);
//	pileup->SaveAs("pileup_MC_RunIISpring15DR74-Asympt25ns.root");
		pileup->SaveAs("nominal.root");
		
}
