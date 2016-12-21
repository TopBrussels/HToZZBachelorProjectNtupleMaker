import os, sys, getopt, argparse,array
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
parser.add_argument('-x','--xsec',help='cross section, in pb', required=False)
parser.add_argument('-l','--lumi',help='luminosity, in pb', required=False)

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
xsec=1
lumi=35000
if args.nMuons :
    nMuo=float(args.nMuons)
    skimmuo=True
if args.nElectrons :
    nEle=float(args.nElectrons)
    skimele=True
if args.tree :
    treename=args.tree
if args.input :
    infile=args.input
if args.output :
    outname=args.output
if args.pTcut :
    minpt=float(args.pTcut)
if args.xsec :
    xsec=float(args.xsec)
if args.lumi :
    lumi=float(args.lumi)
print skimmuo,skimele,nMuo,nEle,treename,infile,outname,minpt


# the analysis structure see TTree/TChain description on root.cern.ch
ch = ROOT.TChain(treename,treename)
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
ch.Add(infile)
# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()

newfile = ROOT.TFile(outname,"recreate")
newtree = ch.CloneTree(0)

bookkeepingtree = ROOT.TChain("startevents","startevents")
bookkeepingtree.Add(infile)

skimweight=array.array("d",[0.0])
skimweight[0]=xsec*lumi/bookkeepingtree.GetEntries()
newbranch = newtree.Branch("skimweight",skimweight,"skimweight/D")

print skimweight

# book some histograms
# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()
lve2=ROOT.TLorentzVector()
# for bookkeeping
ii=0
nevents=ch.GetEntries()
ntimes=100000
# start of loop over events
for iev in ch:
    
    if ii % ntimes ==0 and ii!= 0 :
        print ii, "/", nevents, " fraction skimmed ",100*newtree.GetEntries()/ii,"%"
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
#        print iev.nElectrons
        if iev.nElectrons < nEle:
            continue
#        print iev.nElectrons
        ngoodele=0
        for iele in range(0,iev.nElectrons):
#            print iev.pT_electron[iele]
            if iev.pT_electron[iele] >= minpt :
                ngoodele+=1
        if ngoodele < nEle :
            continue
        newtree.Fill()
    elif skimmuo == True :
#        print iev.nMuons
        if iev.nMuons < nMuo:
            continue
#        print iev.nMuons
        ngoodmuo=0
        for imu in range(0,iev.nMuons):
#            print iev.pT_muon[imu]
            if iev.pT_muon[imu] >= minpt :
                ngoodmuo+=1
        if ngoodmuo < nMuo :
            continue
        newtree.Fill()

# end of loop
newtree.Print()
newtree.AutoSave()
