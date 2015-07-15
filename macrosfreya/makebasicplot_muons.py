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

lumi=11.8
# book some histograms
outfile = ROOT.TFile("output_muons.root","recreate")
outfile.cd()
b_mupt = ROOT.TH1F("b_mupt","muon p_{T}",100,0,500)
b_mueta = ROOT.TH1F("b_mueta","muon #eta",100,-4,4)
b_mud0 = ROOT.TH1F("b_mud0","muon |d_{0}|",200,0,2)
b_muIso = ROOT.TH1F("b_muIso",":muon pfIso:muons",100,0,1)
b_njets = ROOT.TH1F("b_njets","jet mult:jet mult:events",10,0,10)
b_jetpt = ROOT.TH1F("b_jetpt","jet p_{T}",100,0,500)
b_zpeak = ROOT.TH1F("b_zpeak","m(mumu)",100,0,200)


stack_mupt = ROOT.THStack("stack_mupt","muon p_{T}")
#stack_mupt.SetDirectory(outfile)
stack_mueta = ROOT.THStack("stack_mueta","muon #eta")
#stack_mueta.SetDirectory(outfile)
stack_muIso = ROOT.THStack("stack_muIso","muon pfIso")
#stack_muIso.SetDirectory(outfile)
stack_mud0 = ROOT.THStack("stack_mud0","muon d0")
#stack_mud0.SetDirectory(outfile)
stack_munjets = ROOT.THStack("stack_munjets","jet mult:jet mult:events")
#stack_njets.SetDirectory(outfile)
stack_mujetpt = ROOT.THStack("stack_mujetpt","jet p_{T}")
#stack_jetpt.SetDirectory(outfile)
stack_muzpeak = ROOT.THStack("stack_zpeak","m(mumu)")
#stack_zpeak.SetDirectory(outfile)

names=["Wjets","Zjets","ttbar","data"]
colors=[ROOT.kGreen-3,ROOT.kAzure-2,ROOT.kRed+1,ROOT.kBlack]
xsecs=[20508.9,2008.4,800,-1]
filenames=["../../../datafiles/CMSSW_74X_v2-ntuples/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-50ns*.root",
        "../../../datafiles/CMSSW_74X_v2-ntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-50n*.root",
        "../../../datafiles/CMSSW_74X_v2-ntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-50ns*.root",
        "../../../datafiles/CMSSW_74X_v2-ntuples/SingleMuon*.root"]

print xsecs
print filenames

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
    ch.Add(filenames[isam])
    outfile.cd()
    h_mupt=b_mupt.Clone("h_mupt_"+names[isam])
    h_mueta=b_mueta.Clone("h_mueta_"+names[isam])
    h_mud0=b_mud0.Clone("h_mud0_"+names[isam])
    h_muIso=b_muIso.Clone("h_muIso_"+names[isam])
    h_munjets=b_njets.Clone("h_munjets_"+names[isam])
    h_mujetpt=b_jetpt.Clone("h_mujetpt_"+names[isam])
    h_muzpeak=b_zpeak.Clone("h_muzpeak"+names[isam])
    

    
    ii=0
    nevents=ch.GetEntries()
    workxsec=xsecs[isam]
    mcweight=workxsec*lumi
    mcweight/=nevents
    if workxsec == -1:
        mcweight=1  # because then it's data!
    h_mupt.SetLineColor(colors[isam])
    h_mueta.SetLineColor(colors[isam])
    h_mud0.SetLineColor(colors[isam])
    h_muIso.SetLineColor(colors[isam])
    h_munjets.SetLineColor(colors[isam])
    h_mujetpt.SetLineColor(colors[isam])
    h_muzpeak.SetLineColor(colors[isam])
    if workxsec != -1:
        h_mupt.SetFillColor(colors[isam])
        h_mueta.SetFillColor(colors[isam])
        h_mud0.SetFillColor(colors[isam])
        h_muIso.SetFillColor(colors[isam])
        h_munjets.SetFillColor(colors[isam])
        h_mujetpt.SetFillColor(colors[isam])
        h_muzpeak.SetFillColor(colors[isam])



    print "now running on sample ",filenames[isam]," with xsec:",xsecs[isam]," and lumi ",lumi, " and events: ", nevents," giving an MC weight of",mcweight
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
        if iev.nElectrons > 0 :
            continue
        
# loop over muons - fill in lorentz vector and fill some histograms
        for imu in range(0,iev.nMuons) :
            
            lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
            if imu > 0 :
                lve.SetPxPyPzE(iev.pX_muon[imu-1],iev.pY_muon[imu-1],iev.pZ_muon[imu-1],iev.E_muon[imu-1])
                #                print "muon 0: ",lve.Pt()," muon 1:", lvmu.Pt()
                h_muzpeak.Fill((lve+lvmu).M(),mcweight)
            h_mupt.Fill(lvmu.Pt(),mcweight)
            h_mueta.Fill(lvmu.Eta(),mcweight)
            h_mud0.Fill(abs(iev.d0_muon[imu]),mcweight)
            h_muIso.Fill(iev.pfIso_muon[imu],mcweight)



        if iev.nMuons == 1 :
            h_munjets.Fill(iev.nJets,mcweight)
            ntags=0;
            for ijet in range(0, iev.nJets) :
                lvjet.SetPxPyPzE(iev.pX_jet[ijet],iev.pY_jet[ijet],iev.pZ_jet[ijet],iev.E_jet[ijet])
                h_mujetpt.Fill(lvjet.Pt(),mcweight)

                             
    #end of sample loop
    if xsecs[isam] != -1 :
        stack_mujetpt.Add(h_mujetpt)
        stack_munjets.Add(h_munjets)
        stack_mupt.Add(h_mupt)
        stack_mueta.Add(h_mueta)
        stack_muzpeak.Add(h_muzpeak)
    
    outfile.Write()
    print "ran over ", ii, " events"
# end of all loops
print runlist
#
## create canvas
t3=ROOT.TCanvas()
stack_munjets.Draw()
#h_mupt_data.Draw("psame")
#
## create sub-pads and cd() to them, draw some histograms
#t3.Divide(2,2)
#t3.cd(1)
#ROOT.gPad.SetLogy()
#h_mueta.Draw()
#t3.cd(2)
#ROOT.gPad.SetLogy()
#h_mupt.Draw()
#t3.cd(3)
#ROOT.gPad.SetLogy()
#h_mud0.Draw()
#t3.cd(4)
#ROOT.gPad.SetLogy()
#h_njets.Draw()
#t3.Update()
#t3.Print("basic_muonstuff_mujetssel.gif") # TCanvas::Print also can make .pdf or .root/.C plots
#                      
