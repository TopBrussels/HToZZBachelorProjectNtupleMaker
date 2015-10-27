import os, sys
import ROOT
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

# the analysis structure see TTree/TChain description on root.cern.ch
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!

# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()

lumi=1000.0
# book some histograms
outfile = ROOT.TFile("output_electron.root","recreate")
outfile.cd()
b_elept = ROOT.TH1F("b_elept","electron p_{T}",100,0,500)
b_elept.SetXTitle("electron p_{T} [GeV]")
b_eleeta = ROOT.TH1F("b_eleeta","electron #eta",100,-4,4)
b_eleeta.SetXTitle("electron #eta")
b_eled0 = ROOT.TH1F("b_eled0","electron |d_{0}|",200,0,2)
b_eled0.SetXTitle("electron d_{0} [cm]")
b_eleIso = ROOT.TH1F("b_eleIso",":electron pfIso:electrons",100,0,1)
b_eleIso.SetXTitle("electron relative pfIso")
b_elenjets = ROOT.TH1F("b_elenjets","jet elelt:jet elelt:events",10,0,10)
b_elenjets.SetXTitle("jet multiplicity")
b_elejetpt = ROOT.TH1F("b_elejetpt","jet p_{T}",100,0,500)
b_elejetpt.SetXTitle("jet p_{T} [GeV]")
b_elezpeak = ROOT.TH1F("b_elezpeak","m(eleele)",100,0,200)
b_elezpeak.SetXTitle("M(ee) [GeV/]")
b_elenvtx = ROOT.TH1F("b_elenvtx","PV multiplicity",40,0,40)
b_elenvtx.SetXTitle("PV multiplicity")

stack_elept = ROOT.THStack("stack_elept","electron p_{T}")
#stack_elept.SetDirectory(outfile)
stack_eleeta = ROOT.THStack("stack_eleeta","electron #eta")
#stack_eleeta.SetDirectory(outfile)
stack_eleIso = ROOT.THStack("stack_eleIso","electron pfIso")
#stack_eleIso.SetDirectory(outfile)
stack_eled0 = ROOT.THStack("stack_eled0","electron d0")
#stack_eled0.SetDirectory(outfile)
stack_elenjets = ROOT.THStack("stack_elenjets","jet elelt:jet elelt:events")
#stack_njets.SetDirectory(outfile)
stack_elejetpt = ROOT.THStack("stack_elejetpt","jet p_{T}")
#stack_jetpt.SetDirectory(outfile)
stack_elezpeak = ROOT.THStack("stack_elezpeak","m(eleele)")
#stack_zpeak.SetDirectory(outfile)
stack_nvtx = ROOT.THStack("stack_nvtx","PV mult");

names=["Wjets","Zjets","ttbar","tw","atw","data"]
colors=[ROOT.kGreen-3,ROOT.kAzure-2,ROOT.kRed+1,ROOT.kPink,ROOT.kPink,ROOT.kBlack]
xsecs=[20508.9*3,2008.4*3,831.76,35.85,35.85,-1]
filenames=["../WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Ntupler/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_*.root",
           "../DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Ntupler/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_*.root",
           "../TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Ntupler/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_*.root",
           "../ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_Ntupler/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_*.root",
           "../ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_Ntupler/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_*.root",
           "SingleElectron*.root"]

startevcounts=[0,0,0,0,0,0,0,0]
ntupleevcounts=[0,0,0,0,0,0,0,0]


# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()
lvjet=ROOT.TLorentzVector()


runlist=[]
#print "starting, with ", runlist.size() , " runs"

# for bookkeeping

# start of loop over events
for isam in range(len(xsecs)) :
    ch = ROOT.TChain("tree","tree")
#    ch = ROOT.TChain("dileptree","dileptree")
    ch.Add(filenames[isam])
    
    bookkeeping = ROOT.TChain("startevents","startevents")
    bookkeeping.Add(filenames[isam])

#    if isam == 1 :
#        ch.Print()
    outfile.cd()
    h_elept=b_elept.Clone("h_elept_"+names[isam])
    h_eleeta=b_eleeta.Clone("h_eleeta_"+names[isam])
    h_eled0=b_eled0.Clone("h_eled0_"+names[isam])
    h_eleIso=b_eleIso.Clone("h_eleIso_"+names[isam])
    h_elenjets=b_elenjets.Clone("h_elenjets_"+names[isam])
    h_elejetpt=b_elejetpt.Clone("h_elejetpt_"+names[isam])
    h_elezpeak=b_elezpeak.Clone("h_elezpeak_"+names[isam])
    h_elenvtx=b_elenvtx.Clone("h_elenvtx_"+names[isam])
    h_elenopunvtx=b_elenvtx.Clone("h_elenopunvtx_"+names[isam])
    

    
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
    if workxsec == -1:
        mcweight=1  # because then it's data!
    h_elept.SetLineColor(colors[isam])
    h_eleeta.SetLineColor(colors[isam])
    h_eled0.SetLineColor(colors[isam])
    h_eleIso.SetLineColor(colors[isam])
    h_elenjets.SetLineColor(colors[isam])
    h_elejetpt.SetLineColor(colors[isam])
    h_elezpeak.SetLineColor(colors[isam])
    h_elenvtx.SetLineColor(colors[isam])
    h_elenopunvtx.SetLineColor(colors[isam])
    if workxsec != -1:
        h_elept.SetFillColor(colors[isam])
        h_eleeta.SetFillColor(colors[isam])
        h_eled0.SetFillColor(colors[isam])
        h_eleIso.SetFillColor(colors[isam])
        h_elenjets.SetFillColor(colors[isam])
        h_elejetpt.SetFillColor(colors[isam])
        h_elezpeak.SetFillColor(colors[isam])
        h_elenvtx.SetFillColor(colors[isam])
        h_elenopunvtx.SetFillColor(colors[isam])



    print "now running on sample ",filenames[isam]," with xsec:",xsecs[isam]," and lumi ",lumi, " and events: ", nevents," (started with ",neventsstart,"),  giving an MC weight of",mcweight
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
        ngoodjets=0
        if iev.nMuons > 0 :
            continue
        if abs(iev.mc_baseweight - 1.) > 0.00001 :
			if ii < 100 :
				print "MC weights: PU",iev.pu_weight," and generator: ",iev.mc_baseweight
        totalweight=mcweight*iev.pu_weight
        totalweightnopu=mcweight
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
            if iele > 0 :
                if iev.tight_electron[iele-1] ==0 :
                    continue
                lve.SetPtEtaPhiE(iev.pT_electron[iele-1],iev.eta_electron[iele-1],iev.phi_electron[iele-1],iev.E_electron[iele-1])
				
                h_elezpeak.Fill((lve+lvmu).M(),totalweight)
            h_elept.Fill(lvmu.Pt(),totalweight)
            h_eleeta.Fill(lvmu.Eta(),totalweight)
            h_eled0.Fill(abs(iev.d0_electron[iele]),totalweight)
            h_eleIso.Fill(iev.pfIso_electron[iele],totalweight)

#        print "mc weight = ",totalweight    

        if ngoodelectrons == 1:
            h_elenjets.Fill(iev.nJets,totalweight)
            h_elenvtx.Fill(iev.nvtx,totalweight)
            h_elenopunvtx.Fill(iev.nvtx,totalweightnopu)
            ntags=0;
            for ijet in range(0, iev.nJets) :
                lvjet.SetPtEtaPhiE(iev.pT_jet[ijet],iev.eta_jet[ijet],iev.phi_jet[ijet],iev.E_jet[ijet])
                h_elejetpt.Fill(lvjet.Pt(),totalweight)

                             
    #end of sample loop
    if xsecs[isam] != -1 :
        stack_elejetpt.Add(h_elejetpt)
        stack_elenjets.Add(h_elenjets)
        stack_elept.Add(h_elept)
        stack_eleeta.Add(h_eleeta)
        stack_elezpeak.Add(h_elezpeak)
    
    outfile.Write()
    print "ran over ", ii, " events"
# end of all loops
print runlist

print names
print ntupleevcounts
print startevcounts
print xsecs
#
## create canvas
t3=ROOT.TCanvas()
stack_elenjets.Draw("hist")
h_elenjets.Draw("pe0same")
