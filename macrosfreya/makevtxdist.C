{
  TH1F *pileup = new TH1F("pileup","pileup",100,0,100);
  pileup->Sumw2();

  TChain *ch = new TChain("startevents","startevents");
  ch->Add("DoubleEG/*.root");
  ch->Add("DoubleMuon/*.root");
  ch->Add("MuonEG/*.root");
  ch->Add("SingleElectron/*.root");
  ch->Add("SingleMuon/*.root");
  std::cout << "using " << ch->GetEntries() << " events to make pu histogram:" << std::endl;
  ch->Draw("nvtx>>pileup");
    
  double scaler = 1.;
  if(ch->GetEntries()>0)
    scaler /= ch->GetEntries();
  pileup->Scale(scaler);
  pileup->SaveAs("pileup_data_Summer15.root");

}
