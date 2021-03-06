import os, sys
import ROOT
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

# the analysis structure see TTree/TChain description on root.cern.ch
ch = ROOT.TChain("tree","tree")
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
ch.Add("../../../datafiles/CMSSW_74X_v2-ntuples/MuonEG*.root")
# very loud but useful to know what variables are stored in a tree... it prints them all
ch.Print()

# book some histograms
h_mupt = ROOT.TH1F("h_mupt","muon p_{T}",100,0,500)
h_mueta = ROOT.TH1F("h_mueta","muon #eta",100,-4,4)
h_elept = ROOT.TH1F("h_elept","ele p_{T}",100,0,500)
h_eleeta = ROOT.TH1F("h_eleeta","ele #eta",100,-4,4)
h_mud0 = ROOT.TH1F("h_mud0","muon |d_{0}|",200,0,2)
h_eled0 = ROOT.TH1F("h_eled0","ele |d_{0}|",200,0,2)

h_njets = ROOT.TH1F("h_njets","jet mult",10,0,10)

# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()
lvjet=ROOT.TLorentzVector()

# for bookkeeping
ii=0
nevents=ch.GetEntries()

# start of loop over events
for iev in ch:
    if ii % 10000 ==0 :
        print ii, "/", nevents
    ii+=1
# comment out the following lines if you are testing, it stops after a certain number of events
#    if ii==10000 :
#        break


    if iev.nMuons ==0 :
        continue
    if iev.nElectrons ==0 :
        continue
    h_njets.Fill(iev.nJets)
# loop over muons - fill in lorentz vector and fill some histograms
    for imu in range(0,iev.nMuons) :
        lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
        
        h_mupt.Fill(lvmu.Pt())
        h_mueta.Fill(lvmu.Eta())
        h_mud0.Fill(abs(iev.d0_muon[imu]))

# loop over electrons - fill in lorentz vector and fill some histograms
    for iele in range(0,iev.nElectrons) :

        lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])
        h_eleeta.Fill(lve.Eta())
        h_elept.Fill(lve.Pt())
        h_eled0.Fill(abs(iev.d0_electron[iele]))
    if iev.nJets >= 2:
        print  iev.run_num,":",iev.lumi_num,":",iev.evt_num
        print "electron Pt:",lve.Pt()," eta:",lve.Eta()," phi:",lve.Phi()
        print "muon Pt:",lvmu.Pt()," eta:",lvmu.Eta()," phi:",lvmu.Phi()
        print "met: ",iev.missingEt
        for ijet in range(0,iev.nJets) :
            lvjet.SetPxPyPzE(iev.pX_jet[ijet],iev.pY_jet[ijet],iev.pZ_jet[ijet],iev.E_jet[ijet])
            print "jet Pt:",lvjet.Pt()," eta:",lvjet.Eta()," phi:",lvjet.Phi()," btag:",iev.btag_jet[ijet]

# end of loop


# create canvas
t3=ROOT.TCanvas()

# create sub-pads and cd() to them, draw some histograms
t3.Divide(3,2)
t3.cd(1)
h_mueta.Draw()
t3.cd(2)
h_mupt.Draw()
t3.cd(3)
h_eleeta.Draw()
t3.cd(4)
h_elept.Draw()
t3.cd(5)
ROOT.gPad.SetLogy()
h_mud0.Draw()
t3.cd(6)
ROOT.gPad.SetLogy()
h_eled0.Draw()
t3.Update()
t3.Print("basicleptonID.gif") # TCanvas::Print also can make .pdf or .root/.C plots

t4=ROOT.TCanvas()
t4.cd()
h_njets.Draw()
t4.Print("jetmult.gif")