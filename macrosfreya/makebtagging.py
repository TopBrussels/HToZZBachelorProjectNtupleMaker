import os, sys
import ROOT
import time
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
outfile = ROOT.TFile("output_btagging_"+datestring+".root","recreate")
outfile.cd()
b_jetPtEta = ROOT.TH2D("b_jetPtEta","",100,0,500,10,-2.5,2.5)
b_jetPtEta.Sumw2()
b_jetPtEta_HF = b_jetPtEta.Clone("b_jetPtEta_HF")
b_jetPtEta_LF = b_jetPtEta.Clone("b_jetPtEta_LF")
b_jetPtEta_C = b_jetPtEta.Clone("b_jetPtEta_C")
b_jetPtEta_Loose_HF = b_jetPtEta.Clone("b_jetPtEta_Loose_HF")
b_jetPtEta_Loose_LF = b_jetPtEta.Clone("b_jetPtEta_Loose_LF")
b_jetPtEta_Loose_C = b_jetPtEta.Clone("b_jetPtEta_Loose_C")
b_jetPtEta_Medium_HF = b_jetPtEta.Clone("b_jetPtEta_Medium_HF")
b_jetPtEta_Medium_LF = b_jetPtEta.Clone("b_jetPtEta_Medium_LF")
b_jetPtEta_Medium_C = b_jetPtEta.Clone("b_jetPtEta_Medium_C")
b_jetPtEta_Tight_HF = b_jetPtEta.Clone("b_jetPtEta_Tight_HF")
b_jetPtEta_Tight_LF = b_jetPtEta.Clone("b_jetPtEta_Tight_LF")
b_jetPtEta_Tight_C = b_jetPtEta.Clone("b_jetPtEta_Tight_C")
efficiency_Loose_LF = b_jetPtEta.Clone("efficiency_Loose_LF")
efficiency_Medium_LF = b_jetPtEta.Clone("efficiency_Medium_LF")
efficiency_Tight_LF = b_jetPtEta.Clone("efficiency_Tight_LF")
efficiency_Loose_HF = b_jetPtEta.Clone("efficiency_Loose_HF")
efficiency_Medium_HF = b_jetPtEta.Clone("efficiency_Medium_HF")
efficiency_Tight_HF = b_jetPtEta.Clone("efficiency_Tight_HF")
efficiency_Loose_C = b_jetPtEta.Clone("efficiency_Loose_C")
efficiency_Medium_C = b_jetPtEta.Clone("efficiency_Medium_C")
efficiency_Tight_C = b_jetPtEta.Clone("efficiency_Tight_C")



names=["ttbar","tw","atw","tttt"]
colors=[ROOT.kGreen-3,ROOT.kAzure-2,ROOT.kRed+1,ROOT.kPink,ROOT.kPink,ROOT.kGray,ROOT.kBlack]
xsecs=[831.76,35.85,35.85,0.009]
filenames=["ntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia*.root",
           "ntuples/ST_tW_top_5f_inclusiveDecays_13TeV-powheg*.root",
           "ntuples/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-*.root",
           "ntuples/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",
      ]

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
    ch = ROOT.TChain("tree","tree")
#    ch = ROOT.TChain("dileptree","dileptree")
    ch.Add(filenames[isam])
    
    bookkeeping = ROOT.TChain("startevents","startevents")
    bookkeeping.Add(filenames[isam])

#    if isam == 1 :
#        ch.Print()
    outfile.cd()
    

    
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

        if abs(iev.mc_baseweight - 1.) > 0.01 or ii< 100 :
            print "MC weights: PU",iev.pu_weight," and generator: ",iev.mc_baseweight
        totalweight=mcweight*iev.pu_weight*iev.mc_baseweight
        totalweightnopu=mcweight*iev.mc_baseweight
        
        # jet requirement. Needed to get similar kinematics/occupancy as in analysis
        if iev.nJets < 4 :
            continue
        
        ngoodjets=0
        lepweight=1
        ntags=0
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

        ngoodmuons = 0
        for imu in range(0, iev.nMuons):
            ngoodmuons+=1
            lepweight*=iev.sfl_muon[imu]
    
        
        if ngoodelectrons+ngoodmuons != 1 :
            continue
#        print "mc weight = ",totalweight

        for ijet in range(0, iev.nJets) :
            
            totalweight=lepweight*totalweight
            totalweightnopu=lepweight*totalweightnopu
            if iev.flav_jet[ijet] < 4:
                b_jetPtEta_LF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.605 :
                    b_jetPtEta_Loose_LF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.890 :
                    b_jetPtEta_Medium_LF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.970 :
                    b_jetPtEta_Tight_LF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)

            if iev.flav_jet[ijet] == 4:
                b_jetPtEta_C.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.605 :
                    b_jetPtEta_Loose_C.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.890 :
                    b_jetPtEta_Medium_C.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.970 :
                    b_jetPtEta_Tight_C.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
        

            if iev.flav_jet[ijet] == 5:
                b_jetPtEta_HF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.605 :
                    b_jetPtEta_Loose_HF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.890 :
                    b_jetPtEta_Medium_HF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                if iev.btag_jet[ijet] > 0.970 :
                    b_jetPtEta_Tight_HF.Fill(iev.pT_jet[ijet],iev.eta_jet[ijet],totalweight)
                


        storedevents[isam]+=totalweight
    outfile.Write("kOverwrite")
    print "ran over ", ii, " events, selected ",storedevents[isam]
# end of all loops
efficiency_Loose_HF.Divide(b_jetPtEta_Loose_HF,b_jetPtEta_HF,1.,1.,"B")
efficiency_Medium_HF.Divide(b_jetPtEta_Medium_HF,b_jetPtEta_HF,1.,1.,"B")
efficiency_Tight_HF.Divide(b_jetPtEta_Tight_HF,b_jetPtEta_HF,1.,1.,"B")

efficiency_Loose_LF.Divide(b_jetPtEta_Loose_LF,b_jetPtEta_LF,1.,1.,"B")
efficiency_Medium_LF.Divide(b_jetPtEta_Medium_LF,b_jetPtEta_LF,1.,1.,"B")
efficiency_Tight_LF.Divide(b_jetPtEta_Tight_LF,b_jetPtEta_LF,1.,1.,"B")

efficiency_Loose_C.Divide(b_jetPtEta_Loose_C,b_jetPtEta_C,1.,1.,"B")
efficiency_Medium_C.Divide(b_jetPtEta_Medium_C,b_jetPtEta_C,1.,1.,"B")
efficiency_Tight_C.Divide(b_jetPtEta_Tight_C,b_jetPtEta_C,1.,1.,"B")


for ibin in range(1,b_jetPtEta_HF.GetNbinsX()+1) :
    for jbin in range(1,b_jetPtEta_HF.GetNbinsY()+1) :
        print "efficiency for ibin,jbin :  ",ibin,jbin," is L:",efficiency_Loose_HF.GetBinContent(ibin,jbin),"+/-",efficiency_Loose_HF.GetBinError(ibin,jbin), " M:",efficiency_Medium_HF.GetBinContent(ibin,jbin),"+/-",efficiency_Medium_HF.GetBinError(ibin,jbin), " T:", efficiency_Tight_HF.GetBinContent(ibin,jbin),"+/-",efficiency_Tight_HF.GetBinError(ibin,jbin)


#print runlist
for isam in range(len(names)) :
    print "sample ",names[isam]," ",ntupleevcounts[isam]," ",startevcounts[isam]," ",storedevents[isam]," ",xsecs[isam]
#print names
#print ntupleevcounts
#print startevcounts
#print storedevents
#print xsecs
