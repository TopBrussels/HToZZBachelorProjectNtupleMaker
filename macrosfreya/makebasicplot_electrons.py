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
outfile = ROOT.TFile("output_electron_"+datestring+".root","recreate")
outfile.cd()
b_elept = ROOT.TH1F("b_elept","electron p_{T}",100,0,500)
b_elept.SetXTitle("electron p_{T} [GeV]")
b_eleeta = ROOT.TH1F("b_eleeta","electron #eta",100,-4,4)
b_eleeta.SetXTitle("electron #eta")
b_eled0 = ROOT.TH1F("b_eled0","electron |d_{0}|",200,0,2)
b_eled0.SetXTitle("electron d_{0} [cm]")
b_eleDR = ROOT.TH1F("b_eleDR","electron dR wrt closest jet",200,0,2)
b_eleDR.SetXTitle("electron jet #Delta R")
b_eled0zoom = ROOT.TH1F("b_eled0zoom","electron |d_{0}|",200,0,0.02)
b_eled0zoom.SetXTitle("electron d_{0} [cm]")
b_eleIso = ROOT.TH1F("b_eleIso",":electron pfIso:electrons",100,0,1)
b_eleIso.SetXTitle("electron relative pfIso")
b_elenjets = ROOT.TH1F("b_elenjets","jet elelt:jet elelt:events",10,0,10)
b_elenjets.SetXTitle("jet multiplicity")
b_elejetpt = ROOT.TH1F("b_elejetpt","jet p_{T}",100,0,500)
b_elejetpt.SetXTitle("jet p_{T} [GeV]")
b_elezpeak = ROOT.TH1F("b_elezpeak","m(eleele)",100,0,200)
b_elezpeak.SetXTitle("M(ee) [GeV]")
b_elenvtx = ROOT.TH1F("b_elenvtx","PV multiplicity",40,0,40)
b_elenvtx.SetXTitle("PV multiplicity")
b_eleht = ROOT.TH1F("b_eleht","",100,0,1000)
b_eleht.SetXTitle("H_{T} [GeV]")
b_elehtbinned=ROOT.TH1F("b_elehtbinned","",60,0,6000)
b_elehtbinned.SetXTitle("H_{T} [GeV]")


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
           "ntuples/SingleElectron-Run2015*.root"]
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
    h_elept=b_elept.Clone("h_elept_"+names[isam])
    h_eleeta=b_eleeta.Clone("h_eleeta_"+names[isam])
    h_eled0=b_eled0.Clone("h_eled0_"+names[isam])
    h_eleDR=b_eleDR.Clone("h_eleDR_"+names[isam])
    h_eled0zoom=b_eled0zoom.Clone("h_eled0zoom_"+names[isam])
    h_eleIso=b_eleIso.Clone("h_eleIso_"+names[isam])
    h_elenjets=b_elenjets.Clone("h_elenjets_"+names[isam])
    h_elejetpt=b_elejetpt.Clone("h_elejetpt_"+names[isam])
    h_elezpeak=b_elezpeak.Clone("h_elezpeak_"+names[isam])
    h_elenvtx=b_elenvtx.Clone("h_elenvtx_"+names[isam])
    h_elenopunvtx=b_elenvtx.Clone("h_elenopunvtx_"+names[isam])
    h_eleht=b_eleht.Clone("h_eleht_"+names[isam])
    h_elehtbinned=b_elehtbinned.Clone("h_elehtbinned_"+names[isam])

    

    
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
    h_elept.SetLineColor(colors[isam])
    h_eleeta.SetLineColor(colors[isam])
    h_eleDR.SetLineColor(colors[isam])
    h_eled0.SetLineColor(colors[isam])
    h_eled0zoom.SetLineColor(colors[isam])
    h_eleIso.SetLineColor(colors[isam])
    h_elenjets.SetLineColor(colors[isam])
    h_elejetpt.SetLineColor(colors[isam])
    h_elezpeak.SetLineColor(colors[isam])
    h_elenvtx.SetLineColor(colors[isam])
    h_elenopunvtx.SetLineColor(colors[isam])
    h_eleht.SetLineColor(colors[isam])
    h_elehtbinned.SetLineColor(colors[isam])

    if workxsec != -1:
        h_elept.SetFillColor(colors[isam])
        h_eleeta.SetFillColor(colors[isam])
        h_eled0.SetFillColor(colors[isam])
        h_eleDR.SetFillColor(colors[isam])
        h_eled0zoom.SetFillColor(colors[isam])
        h_eleIso.SetFillColor(colors[isam])
        h_elenjets.SetFillColor(colors[isam])
        h_elejetpt.SetFillColor(colors[isam])
        h_elezpeak.SetFillColor(colors[isam])
        h_elenvtx.SetFillColor(colors[isam])
        h_elenopunvtx.SetFillColor(colors[isam])
        h_eleht.SetFillColor(colors[isam])
        h_elehtbinned.SetFillColor(colors[isam])



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
        if iev.nMuons > 0 :
            continue
        if iev.nJets < 4 :
            continue
        ntags=0

        ht=0

        #ad-hoc calculation of btag sfs:
        arraywithsfs=[]
        
        for ijet in range(0, iev.nJets) :
            lvjet.SetPtEtaPhiE(iev.pT_jet[ijet],iev.eta_jet[ijet],iev.phi_jet[ijet],iev.E_jet[ijet])
            h_elejetpt.Fill(lvjet.Pt(),totalweight)
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

        ngoodelectrons = 0
        if workxsec == -1 :
            #reset weight for data:
            totalweight=1 
# loop over electrons - fill in lorentz vector and fill some histograms



        for iele in range(0,iev.nElectrons) :
            if iev.tight_electron[iele] == 0 :
                continue
            lvmu.SetPtEtaPhiE(iev.pT_electron[iele],iev.eta_electron[iele],iev.phi_electron[iele],iev.E_electron[iele])
            ngoodelectrons+=1
            lepweight*=iev.sfl_electron[iele];
            if iele > 0 :
                if iev.tight_electron[iele-1] ==0 :
                    continue
                lve.SetPtEtaPhiE(iev.pT_electron[iele-1],iev.eta_electron[iele-1],iev.phi_electron[iele-1],iev.E_electron[iele-1])
                if workxsec != -1 :
                    h_elezpeak.Fill((lve+lvmu).M(),totalweight*lepweight*iev.sfl_electron[iele-1])
                else :
                    h_elezpeak.Fill((lve+lvmu).M(),1)
            if workxsec == -1 :
                h_elept.Fill(lvmu.Pt(),1)
                h_eleeta.Fill(lvmu.Eta(),1)
                h_eled0.Fill(abs(iev.d0_electron[iele]),1)
                h_eled0zoom.Fill(abs(iev.d0_electron[iele]),1)
                h_eleIso.Fill(iev.pfIso_electron[iele],1)
                h_eleDR.Fill(iev.drJet_electron[iele],1)

            else :
                h_elept.Fill(lvmu.Pt(),totalweight*lepweight)
                h_eleeta.Fill(lvmu.Eta(),totalweight*lepweight)
                h_eled0.Fill(abs(iev.d0_electron[iele]),totalweight*lepweight)
                h_eled0zoom.Fill(abs(iev.d0_electron[iele]),totalweight*lepweight)
                h_eleIso.Fill(iev.pfIso_electron[iele],totalweight*lepweight)
                h_eleDR.Fill(iev.drJet_electron[iele],totalweight*lepweight)

        totalweight*=lepweight
        totalweightnopu*=lepweight
#        print "mc weight = ",totalweight
        if ngoodelectrons == 1:
            if workxsec == -1:
                totalweight = 1
                totalweightnopu = 1
            h_elenjets.Fill(iev.nJets,totalweight)
            h_elenvtx.Fill(iev.nvtx,totalweight)
            h_elenopunvtx.Fill(iev.nvtx,totalweightnopu)
            h_eleht.Fill(ht,totalweight)
            htfixed=ht
            if htfixed>1000 :
                htfixed=1000.
            if iev.nJets > 7:
                htfixed+=1000.
            if iev.nJets > 8:
                htfixed+=1000.
            if ntags > 3 :
                htfixed+=3000.
            h_elehtbinned.Fill(htfixed,totalweight)
    storedevents[isam]=h_elenjets.GetSum()
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
