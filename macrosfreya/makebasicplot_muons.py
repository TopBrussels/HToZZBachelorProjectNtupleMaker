import os, sys
import ROOT
import time
import array

import FreyaLib as fl
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

# the analysis structure see TTree/TChain description on root.cern.ch
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!

# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()
lumi=1000.0

datestring=time.strftime("%Y-%m-%d")
print datestring
# book some histograms
outfile = ROOT.TFile("output_muon_"+datestring+".root","recreate")
outfile.cd()
b_muopt = ROOT.TH1F("b_muopt","muon p_{T}",100,0,500)
b_muopt.SetXTitle("muon p_{T} [GeV]")
b_muoeta = ROOT.TH1F("b_muoeta","muon #eta",100,-4,4)
b_muoeta.SetXTitle("muon #eta")
b_muod0 = ROOT.TH1F("b_muod0","muon |d_{0}|",200,0,2)
b_muod0.SetXTitle("muon d_{0} [cm]")
b_muoDR = ROOT.TH1F("b_muoDR","muon dR wrt closest jet",200,0,2)
b_muoDR.SetXTitle("muon jet #Delta R")
b_muod0zoom = ROOT.TH1F("b_muod0zoom","muon |d_{0}|",200,0,0.02)
b_muod0zoom.SetXTitle("muon d_{0} [cm]")
b_muoIso = ROOT.TH1F("b_muoIso",":muon pfIso:muons",100,0,1)
b_muoIso.SetXTitle("muon relative pfIso")
b_muonjets = ROOT.TH1F("b_muonjets","jet muolt:jet muolt:events",10,0,10)
b_muonjets.SetXTitle("jet multiplicity")
b_muojetpt = ROOT.TH1F("b_muojetpt","jet p_{T}",100,0,500)
b_muojetpt.SetXTitle("jet p_{T} [GeV]")
b_muozpeak = ROOT.TH1F("b_muozpeak","m(muomuo)",100,0,200)
b_muozpeak.SetXTitle("M(ee) [GeV]")
b_muonvtx = ROOT.TH1F("b_muonvtx","PV multiplicity",40,0,40)
b_muonvtx.SetXTitle("PV multiplicity")
b_muoht = ROOT.TH1F("b_muoht","",100,0,1000)
b_muoht.SetXTitle("H_{T} [GeV]")
b_muohtbinned=ROOT.TH1F("b_muohtbinned","",60,0,6000)
b_muohtbinned.SetXTitle("H_{T} [GeV]")


skipsample=[0,0,0,0,0,0,0]
names=["Wjets","Zjets","ttbar","tw","atw","tttt","data"]
colors=[ROOT.kGreen-3,ROOT.kAzure-2,ROOT.kRed+1,ROOT.kPink,ROOT.kPink,ROOT.kGray,ROOT.kBlack]
xsecs=[20508.9*3,2008.4*3,831.76,35.85,35.85,0.009,-1]
filenames=["ntuples/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8*.root",
           "ntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-*.root",
           "ntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia*.root",
           "ntuples/ST_tW_top_5f_inclusiveDecays_13TeV-powheg*.root",
           "ntuples/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-*.root",
           "ntuples/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",
           "ntuples/SingleMuon-Run2015*.root"]
sfnegeventratios=[16518190.0997,9052671,11339232,995600,1000000,415094.011903,-1,-1,-1]
startevcounts=[0.,0.,0.,0.,0.,0.,0.,0.,0.]
ntupleevcounts=[0.,0,0.,0.,0.,0.,0.,0.,0.,0.]
storedevents=[0.,0.,0.,0.,0.,0.,0.,0.,0.]


#print runlist
for isam in range(len(names)) :
    print "sample ",names[isam]," ",ntupleevcounts[isam]," ",startevcounts[isam]," ",storedevents[isam]," ",xsecs[isam]
#print names
#print ntupleevcounts


# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()
lvjet=ROOT.TLorentzVector()


runlist=[]
#print "starting, with ", runlist.size() , " runs"

# for bookkeeping

# start of loop over events
for isam in range(len(xsecs)) :
    if skipsample[isam]==1 :
        continue
    ch = ROOT.TChain("tree","tree")
#    ch = ROOT.TChain("dileptree","dileptree")
    ch.Add(filenames[isam])
#    ch.Print()
    
    bookkeeping = ROOT.TChain("startevents","startevents")
    bookkeeping.Add(filenames[isam])

#    if isam == 1 :
#        ch.Print()
    outfile.cd()
    h_muopt=b_muopt.Clone("h_muopt_"+names[isam])
    h_muoeta=b_muoeta.Clone("h_muoeta_"+names[isam])
    h_muod0=b_muod0.Clone("h_muod0_"+names[isam])
    h_muoDR=b_muoDR.Clone("h_muoDR_"+names[isam])
    h_muod0zoom=b_muod0zoom.Clone("h_muod0zoom_"+names[isam])
    h_muoIso=b_muoIso.Clone("h_muoIso_"+names[isam])
    h_muonjets=b_muonjets.Clone("h_muonjets_"+names[isam])
    h_muojetpt=b_muojetpt.Clone("h_muojetpt_"+names[isam])
    h_muozpeak=b_muozpeak.Clone("h_muozpeak_"+names[isam])
    h_muonvtx=b_muonvtx.Clone("h_muonvtx_"+names[isam])
    h_muonopunvtx=b_muonvtx.Clone("h_muonopunvtx_"+names[isam])
    h_muoht=b_muoht.Clone("h_muoht_"+names[isam])
    h_muohtbinned=b_muohtbinned.Clone("h_muohtbinned_"+names[isam])

    

    
    ii=0
    nevents=ch.GetEntries()
    if nevents==0 :
		continue
    neventsstart = bookkeeping.GetEntries()
    startevcounts[isam]=neventsstart
    ntupleevcounts[isam]=nevents
    workxsec=xsecs[isam]
    mcweight=float(workxsec*lumi)
    mcweight/=float(neventsstart)
    lepweight=1.
    if workxsec == -1:
        mcweight=1  # because then it's data!
    h_muopt.SetLineColor(colors[isam])
    h_muoeta.SetLineColor(colors[isam])
    h_muoDR.SetLineColor(colors[isam])
    h_muod0.SetLineColor(colors[isam])
    h_muod0zoom.SetLineColor(colors[isam])
    h_muoIso.SetLineColor(colors[isam])
    h_muonjets.SetLineColor(colors[isam])
    h_muojetpt.SetLineColor(colors[isam])
    h_muozpeak.SetLineColor(colors[isam])
    h_muonvtx.SetLineColor(colors[isam])
    h_muonopunvtx.SetLineColor(colors[isam])
    h_muoht.SetLineColor(colors[isam])
    h_muohtbinned.SetLineColor(colors[isam])

    if workxsec != -1:
        h_muopt.SetFillColor(colors[isam])
        h_muoeta.SetFillColor(colors[isam])
        h_muod0.SetFillColor(colors[isam])
        h_muoDR.SetFillColor(colors[isam])
        h_muod0zoom.SetFillColor(colors[isam])
        h_muoIso.SetFillColor(colors[isam])
        h_muonjets.SetFillColor(colors[isam])
        h_muojetpt.SetFillColor(colors[isam])
        h_muozpeak.SetFillColor(colors[isam])
        h_muonvtx.SetFillColor(colors[isam])
        h_muonopunvtx.SetFillColor(colors[isam])
        h_muoht.SetFillColor(colors[isam])
        h_muohtbinned.SetFillColor(colors[isam])



    mcweightssf=0.0;
    if sfnegeventratios[isam] < 0:
        if xsecs[isam] > 0 :
            for iev in bookkeeping:
                mcweightssf+=iev.mc_baseweight
        else :
            mcweightssf=bookkeeping.GetEntries()
    else :
        mcweightssf =sfnegeventratios[isam]

    print "now running on sample ",filenames[isam]," with xsec:",xsecs[isam]," and lumi ",lumi, " and events: ", nevents," (started with ",neventsstart,"),  giving an MC weight of",mcweight," baseweight: ",mcweightssf
    mcweight*=mcweightssf/neventsstart

    for iev in ch:
        if ii % 10000 ==0 :
            print ii, "/", nevents
        ii+=1
# comment out the following lines if you are testing, it stops after a certain number of events
#    if ii==10000 :
#        break
        if workxsec == -1 :
            if runlist.count(iev.run_num) ==0 :
                runlist.append(iev.run_num);

        if abs(iev.mc_baseweight - 1.) > 0.00001 :
            if ii < 100 :
                print "MC weights: PU",iev.pu_weight," and generator: ",iev.mc_baseweight
        totalweight=mcweight*iev.pu_weight*iev.mc_baseweight
        totalweightnopu=mcweight*iev.mc_baseweight
        
        ngoodjets=0
        lepweight=1
        if iev.nElectrons > 0 :
            continue
        if iev.nJets < 4 :
            continue
        ntags=0

        ht=0

        #ad-hoc calculation of btag sfs:
        arraywithsfs=[]
        
        for ijet in range(0, iev.nJets) :
            lvjet.SetPtEtaPhiE(iev.pT_jet[ijet],iev.eta_jet[ijet],iev.phi_jet[ijet],iev.E_jet[ijet])
            h_muojetpt.Fill(lvjet.Pt(),totalweight)
            ht+=lvjet.Pt()
            arraywithsfs.append(max(0,min(iev.effb_jet[ijet]*iev.sfb_jet[ijet],1.0)))
            #            print "btag sfs:",arraywithoutsfs[ijet],arraywithsfs[ijet]
            if iev.btag_jet[ijet] > 0.890 :
                ntags+=1


        if workxsec ==-1:
            if ntags < 2 :
                continue

        #        totalweight*=iev.mc_btgsfweight2[0]
#        totalweightnopu*=iev.mc_btgsfweight2[0]

        if workxsec == -1:
            totalweight = 1
            totalweightnopu = 1

        
        if workxsec != -1:
            bsfin=1.-fl.p0(arraywithsfs)-fl.p1(arraywithsfs)
            totalweight*=bsfin
            totalweightnopu*=bsfin

        ngoodmuons = 0
        if workxsec == -1 :
            #reset weight for data:
            totalweight=1 
# loop over muons - fill in lorentz vector and fill some histograms



        for imuo in range(0,iev.nMuons) :
            lvmu.SetPtEtaPhiE(iev.pT_muon[imuo],iev.eta_muon[imuo],iev.phi_muon[imuo],iev.E_muon[imuo])
            ngoodmuons+=1
            lepweight*=iev.sfl_muon[imuo];
            if imuo > 0 :
                lve.SetPtEtaPhiE(iev.pT_muon[imuo-1],iev.eta_muon[imuo-1],iev.phi_muon[imuo-1],iev.E_muon[imuo-1])
                if workxsec != -1 :
                    h_muozpeak.Fill((lve+lvmu).M(),totalweight*lepweight*iev.sfl_muon[imuo-1])
                else :
                    h_muozpeak.Fill((lve+lvmu).M(),1)
            if workxsec == -1 :
                h_muopt.Fill(lvmu.Pt(),1)
                h_muoeta.Fill(lvmu.Eta(),1)
                h_muod0.Fill(abs(iev.d0_muon[imuo]),1)
                h_muod0zoom.Fill(abs(iev.d0_muon[imuo]),1)
                h_muoIso.Fill(iev.pfIso_muon[imuo],1)
                h_muoDR.Fill(iev.drJet_muon[imuo],1)

            else :
                h_muopt.Fill(lvmu.Pt(),totalweight*lepweight)
                h_muoeta.Fill(lvmu.Eta(),totalweight*lepweight)
                h_muod0.Fill(abs(iev.d0_muon[imuo]),totalweight*lepweight)
                h_muod0zoom.Fill(abs(iev.d0_muon[imuo]),totalweight*lepweight)
                h_muoIso.Fill(iev.pfIso_muon[imuo],totalweight*lepweight)
                h_muoDR.Fill(iev.drJet_muon[imuo],totalweight*lepweight)

        totalweight*=lepweight
        totalweightnopu*=lepweight
#        print "mc weight = ",totalweight
        if ngoodmuons == 1:
            if workxsec == -1:
                totalweight = 1
                totalweightnopu = 1
            h_muonjets.Fill(iev.nJets,totalweight)
            h_muonvtx.Fill(iev.nvtx,totalweight)
            h_muonopunvtx.Fill(iev.nvtx,totalweightnopu)
            h_muoht.Fill(ht,totalweight)
            htfixed=ht
            if htfixed>1000 :
                htfixed=1000.
            if iev.nJets > 7:
                htfixed+=1000.
            if iev.nJets > 8:
                htfixed+=1000.
            if ntags > 3 :
                htfixed+=3000.
            h_muohtbinned.Fill(htfixed,totalweight)
    storedevents[isam]=h_muonjets.GetSum()
    outfile.Write()
    print "ran over ", ii, " events, selected ",storedevents[isam]
# end of all loops

#print runlist
for isam in range(len(names)) :
    print "sample ",names[isam]," ",ntupleevcounts[isam]," ",startevcounts[isam]," ",storedevents[isam]," ",xsecs[isam]
#print names
#print ntupleevcounts
#print startevcounts
#print storedevents
#print xsecs
