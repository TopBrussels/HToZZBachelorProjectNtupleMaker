import os, sys
import ROOT as rt
import CMS_lumi, tdrstyle
import array


#set the tdr style
tdrstyle.setTDRStyle()

CMS_lumi.lumi_13TeV="L #approx 553 pb^{-1}"
luminosity= 533
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Super-Preliminary Work in Progress"


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

canvas = rt.TCanvas("c2","c2",50,50,W,H)
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

#file=rt.TFile("output_electrons.root","read")
file=rt.TFile("output_muons.root","read")
#file=rt.TFile("output_displaced.root","read")
file.ls()
#file.ls()
#listofnames=["h_eled0_","h_mud0_"]
#listofnames=["h_elenjets_","h_elezpeak_","h_elenvtx_","h_elept_","h_eleeta_","h_eleIso_"]
listofnames=["h_elenvtx_"]

for histname in listofnames :
    print "now making plots: ",histname
    canv =  rt.TCanvas("c2","c2",50,50,W,H)
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
    hist_tt=file.Get(histname+"ttbar")
    #    hist_tt.DrawClone();
    hist_z=file.Get(histname+"Zjets")
    hist_z.Scale(luminosity)
    #    hist_z.DrawClone("same")
    hist_w=file.Get(histname+"Wjets")
    hist_w.Scale(luminosity)
    #    hist_w.DrawClone("same")
    hist_atw=file.Get(histname+"atw")
    hist_atw.Scale(luminosity)
    #    hist_w.DrawClone("same")
    hist_tw=file.Get(histname+"tw")
    hist_tw.Scale(luminosity)
    #    hist_w.DrawClone("same")
    hist_data=file.Get(histname+"data")
    print "number of data events in histogram ",hist_data.GetName()," is ",hist_data.GetIntegral()
    stack = rt.THStack("stack"+histname,"")
    stack.Add(hist_w,"f")
    stack.Add(hist_z,"f")
    stack.Add(hist_tw,"f")
    stack.Add(hist_atw,"f")
    stack.Add(hist_tt,"f")


    canv.cd()
    hist_data.SetYTitle("events")
    hist_data.Draw("ex0")
    
    stack.Draw("histsame")
    leg = rt.TLegend(0.7,0.65,0.95,0.95)
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
    hist_data.Draw("esamex0")
    leg.Draw("same")
    hist_data.SetMinimum(0.01)
    hist_data.SetMaximum(5.*hist_data.GetMaximum())
    rt.gPad.SetLogy()
    canv.cd()
    CMS_lumi.CMS_lumi(canv, 4, iPos)
    canv.Update()
    canv.RedrawAxis()
    frame = canv.GetFrame()
    frame.Draw()
    canv.Print(histname+"stackplot.png")
    canv.Print(histname+"stackplot.C")
raw_input("Type Entry to end")
