import os, sys
import ROOT as rt
import CMS_lumi, tdrstyle
import FreyaLib as fl
import array
import time


datetime=time.strftime("%Y-%m-%d-%H-%M")
print datetime

#set the tdr style
tdrstyle.setTDRStyle()

luminosity= (16344578.402/1000000000.)+(887308177.909/1000000000.)+(1558092799.428/1000000000.) # electrons
#luminosity= (16344578.402/1000000000.)+(886373910.964/1000000000.)#+(1557667314.725/1000000000.) # muons


CMS_lumi.lumi_13TeV="L #approx "+str(luminosity)+" fb^{-1} (SILVER)"
print "luminosity is: ",luminosity
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Work in Progress, e+4jets+2CSVM"
#CMS_lumi.extraText = "Work in Progress, #mu+4jets+2CSVM"


iPos = 10
if( iPos==0 ): CMS_lumi.relPosX = 0.12

H_ref = 600;
W_ref = 800;
W = W_ref
H  = H_ref


T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

#canvas = rt.TCanvas("c2","c2",50,50,W,H)
canvas = fl.makeDivCan()
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)


file=rt.TFile("output_electron_2015-12-04.root","read")
#file=rt.TFile("output_muon_2015-12-04.root","read")
#file=rt.TFile("output_displaced.root","read")
file.ls()
#file.ls()
#listofnames=["h_eled0_","h_mud0_"]
#listofnames=["h_muonvtx_","h_muonopunvtx_","h_muozpeak_","h_muonjets_","h_muopt_","h_muoIso_","h_muod0zoom_","h_muoht_","h_muohtbinned_","h_muoDR_"]
listofnames=["h_elenvtx_","h_elenopunvtx_","h_elezpeak_","h_elenjets_","h_elept_","h_eleIso_","h_eled0zoom_","h_eleht_","h_elehtbinned_","h_eleDR_"]
dirname="plots_"+datetime
rt.gSystem.mkdir(dirname)
for histname in listofnames :
    print "now making plots for: ",histname
    canv = fl.makeDivCan()
    canv.SetFillColor(0)
    canv.SetBorderMode(0)
    canv.SetFrameFillStyle(0)
    canv.SetFrameBorderMode(0)
    canv.SetLeftMargin( L/W )
    canv.SetRightMargin( R/W )
    canv.SetTopMargin( T/H )
    canv.SetBottomMargin( B/H )
    canv.SetTickx(0)
    canv.SetTicky(0)
    print "retrieving ",histname,"ttbar"

    hist_tt=file.Get(histname+"ttbar")
    print "number of entries: ",hist_tt.GetSum()
    hist_tt.Scale(luminosity)
    print "retrieving ",histname,"Zjets"
    hist_z=file.Get(histname+"Zjets")
    print "number of entries: ",hist_z.GetSum()
    hist_z.Scale(luminosity)
    print "retrieving ",histname,"Wjets"
    hist_w=file.Get(histname+"Wjets")
    print "number of entries: ",hist_w.GetSum()
    hist_w.Scale(luminosity)
    print "retrieving ",histname,"atw"
    hist_atw=file.Get(histname+"atw")
    print "number of entries: ",hist_atw.GetSum()
    hist_atw.Scale(luminosity)
    print "retrieving ",histname,"tw"
    hist_tw=file.Get(histname+"tw")
    print "number of entries: ",hist_tw.GetSum()
    hist_tw.Scale(luminosity)
    print "retrieving ",histname,"tttt"
    hist_tttt=file.Get(histname+"tttt")
    print "number of entries: ",hist_tttt.GetSum()
    hist_tttt.Scale(luminosity*100)
    
    print "retrieving ",histname,"data"
    hist_data=file.Get(histname+"data")
    print "number of data events in histogram ",hist_data.GetName()," is ",hist_data.GetSum()
    stack = rt.THStack("stack"+histname,"")
    stack.Add(hist_w,"f")
    stack.Add(hist_z,"f")
    stack.Add(hist_tw,"f")
    stack.Add(hist_atw,"f")
    stack.Add(hist_tt,"f")


    canv.cd(1)
    hist_data.SetYTitle("events")
    hist_data.Draw("ex0")
    hist_tttt.Draw("histsame")
    
    stack.Draw("histsame")
    leg = rt.TLegend(0.7,0.65,0.95,0.90)
    leg.SetFillStyle(0)
    titletext="data :"+"{0:.2f}".format(hist_data.GetSum())
    leg.AddEntry(hist_data,titletext,"p")
    titletext="ttbar :"+"{0:.2f}".format(hist_tt.GetSum())
    leg.AddEntry(hist_tt,titletext,"f")
    titletext="tW :"+"{0:.2f}".format(hist_tw.GetSum()+hist_atw.GetSum())
    leg.AddEntry(hist_tw,titletext,"f")
    titletext="Z+jets :"+"{0:.2f}".format(hist_z.GetSum())
    leg.AddEntry(hist_z,titletext,"f")
    titletext="W+jets :"+"{0:.2f}".format(hist_w.GetSum())
    leg.AddEntry(hist_w,titletext,"f")
    titletext="tttt x100 :"+"{0:.2f}".format(hist_tttt.GetSum())
    leg.AddEntry(hist_tttt,titletext,"l")
    hist_tttt.SetFillStyle(0)
    hist_tttt.SetLineWidth(2*hist_tttt.GetLineWidth())
    hist_tttt.Draw("histsame")
    hist_data.Draw("esamex0")
    leg.Draw("same")
    hist_data.SetMinimum(0.01)
    hist_data.SetMaximum(20*hist_data.GetMaximum())
    rt.gPad.SetLogy()
    canv.cd(1)
    CMS_lumi.CMS_lumi(rt.gPad, 4, iPos)
    canv.Update()
    canv.RedrawAxis()
    frame = canv.GetFrame()
    frame.Draw()
    canv.cd(2)
    rt.gPad.SetGridy()
    ratioplotH = fl.makeRatioDiff(hist_data,stack)
    ratioplotH.SetMaximum(1.)
    ratioplotH.SetMinimum(-1.)
    ratioplotH.Draw()
    line = rt.TLine(ratioplotH.GetBinLowEdge(1),0,ratioplotH.GetBinLowEdge(ratioplotH.GetNbinsX()+1),0)
    line.SetLineColor(rt.kRed)
    line.SetLineStyle(1)
    line.Draw("same")
    canv.cd()
    canv.Print(dirname+"/"+histname+"stackplot"+datetime+".png")
    canv.Print(dirname+"/"+histname+"stackplot"+datetime+".pdf")
    canv.Print(dirname+"/"+histname+"stackplot"+datetime+".C")
raw_input("Type Entry to end")
