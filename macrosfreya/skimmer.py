import os, sys, getopt, argparse
import ROOT
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#
parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i','--input', help='Input file name',required=True)
parser.add_argument('-o','--output', help='Output file name',required=True)
parser.add_argument('-t','--tree', help='Name of tree',required=False)
parser.add_argument('-nMu','--nMuons', help='number of muons cut',required=False)
parser.add_argument('-nEl','--nElectrons', help='number of electrons cut',required=False)
parser.add_argument('-pt','--pTcut',help='min pT cut', required=False)
args = parser.parse_args()

print args

minpt=0
skimele=False
skimmuo=False
nMuo=0
nEle=0
treename="displtree"
outname="dummyfilename.root"
infile="../CMSSW_80X_v1-ntuples/*DoubleEG*.root"

if args.nMuons :
    nMuo=args.nMuons
    skimmuo=True
if args.nElectrons :
    nEle=args.nElectrons
    skimele=True
if args.tree :
    treename=args.tree
if args.input :
    infile=args.input
if args.output :
    outname=args.output
if args.pTcut :
    minpt=args.pTcut
print nMuo,nEle,treename,infile,outname


# the analysis structure see TTree/TChain description on root.cern.ch
ch = ROOT.TChain(treename,treename)
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
ch.Add(infile)
# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()

newfile = ROOT.TFile(outname,"recreate")
newtree = ch.CloneTree(0)

# book some histograms
# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()
lve2=ROOT.TLorentzVector()
# for bookkeeping
ii=0
nevents=ch.GetEntries()

# start of loop over events
for iev in ch:
    if ii % 100000 ==0 :
        print ii, "/", nevents
    ii+=1
# comment out the following lines if you are testing, it stops after a certain number of events
#    if ii==10000 :
#        break


# loop over muons - fill in lorentz vector and fill some histograms
#    for imu in range(0,iev.nMuons) :
#        lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
#        
#        h_mupt.Fill(lvmu.Pt())
#        h_mueta.Fill(lvmu.Eta())
#        h_mud0.Fill(abs(iev.d0_muon[imu]))

# loop over electrons - fill in lorentz vector and fill some histograms
    if skimele == True :
        if iev.nElectrons< nEle:
            continue
        ngoodele=0
        for ii in range(0,iev.nElectrons):
            if iev.pT_electron[ii] >= ptcut :
                ngoodele+=1
        if ngoodele < nEle :
            continue
        newtree.Fill()
    elif skimmuo == True :
        if iev.nMuons < nMuo:
            continue
        ngoodmuo=0
        for ii in range(0,iev.nMuons):
            if iev.pT_muon[ii] >= ptcut :
                ngoodmuo+=1
        if ngoodmuo < nMuo :
            continue
        newtree.Fill()

# end of loop
newtree.Print()
newtree.AutoSave()
