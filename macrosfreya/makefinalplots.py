import os, sys
import ROOT as rt
import CMS_lumi, tdrstyle
import array


#set the tdr style
tdrstyle.setTDRStyle()

CMS_lumi.lumi_13TeV="L #approx 11.8 pb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Super-Preliminary"


iPos = 12
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

file=rt.TFile("output_electrons.root","read")
#file=rt.TFile("output_muons.root","read")
file.ls()
#file.ls()
listofnames=["h_elenjets_","h_elezpeak","h_eled0_"]
#listofnames=["h_munjets_","h_muzpeak","h_mud0_"]

for histname in listofnames :
    print "now making plots: ",histname
    canv =  rt.TCanvas()
    hist_tt=file.Get(histname+"ttbar")
    #    hist_tt.DrawClone();
    hist_z=file.Get(histname+"Zjets")
    #    hist_z.DrawClone("same")
    hist_w=file.Get(histname+"Wjets")
    #    hist_w.DrawClone("same")
    hist_data=file.Get(histname+"data")
    stack = rt.THStack("stack"+histname,"")
    stack.Add(hist_w,"f")
    stack.Add(hist_z,"f")
    stack.Add(hist_tt,"f")


    hist_data.Draw("ex0")
    
    stack.Draw("histsame")
    hist_data.Draw("esamex0")
    hist_data.SetMinimum(0.01)
    rt.gPad.SetLogy()
    canv.cd()
    CMS_lumi.CMS_lumi(canv, 4, iPos)
    canv.Update()

    canv.RedrawAxis()
    frame = canv.GetFrame()
    frame.Draw()
raw_input("Type Entry to end")
