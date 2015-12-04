# lib keeping track of useful pyroot functions

import ROOT

def p0(values):
    retval=1.
    for ival in range(len(values)) :
        retval*=(1.-values[ival])

    return retval

def p1(values):
    retval=0.
    for ival in range(len(values)) :
        #        print values[ival]
        thisweight=values[ival]
        for jval in range(len(values)) :
            if ival == jval :
                continue
        #            print "(1-",values[jval]
            thisweight*=(1.-values[jval])
        retval+=thisweight
#        print thisweight,retval
    return retval

def makeRatio(hist1,hist2,ymax=False,ymin=False,norm=False):
    """returns the ratio plot hist2/hist1
    if one of the histograms is a stack put it in as argument 2!"""

    if norm:
        print 'scaling!'
        try:
            print 'scale 1: ',1/hist1.Integral()
            print 'scale 2: ',1/hist2.Integral()
            hist1.Scale(1/hist1.Integral())
            hist2.Scale(1/hist2.Integral())
        except(ZeroDivisionError):
            pass
    retH = hist1.Clone()
    try:
        retH.Divide(hist2)
    except(TypeError):
        #this is the error you get if hist2 is a stack
        hList = hist2.GetHists()
        sumHist = hist1.Clone("sumHist")
        sumHist.Reset()
        for h in hList:
            sumHist.Add(h)
        retH.Divide(sumHist)
    except(AttributeError):
        #this is the error you get if hist1 is a stack
        print "Did you use a stack as argument 1? please use stack as argument 2!"
        raise AttributeError
    if ymax or ymin:
        retH.GetYaxis().SetRangeUser(0.5,1.5)
        retH.SetLineColor(hist2.GetLineColor())

    retH.GetAxisX().SetLabelSize(hist1.GetHistogram().GetAxisX().GetLabelSize())
    retH.GetAxisX().SetLabelOffset(hist1.GetHistogram().GetAxisX().GetLabelOffset())
    retH.GetAxisX().SetTitleSize(hist1.GetHistogram().GetAxisX().GetTitleSize())
    retH.GetAxisX().GetTitleOffset(hist1.GetHistogram().GetAxisX().GetTitleOffset())
    ROOT.SetOwnership(retH,0)
    return retH


def makeRatioDiff(hist1,hist2,ymax=False,ymin=False,norm=False):
    """returns the ratio plot hist2/hist1
        if one of the histograms is a stack put it in as argument 2!"""
    
    if norm:
        print 'scaling!'
        try:
            print 'scale 1: ',1/hist1.Integral()
            print 'scale 2: ',1/hist2.Integral()
            hist1.Scale(1/hist1.Integral())
            hist2.Scale(1/hist2.Integral())
        except(ZeroDivisionError):
            pass
    retH = hist1.Clone()
    try:
        retH.Add(hist2,-1.)
        retH.Divide(hist2)
    except(TypeError):
        #this is the error you get if hist2 is a stack
        hList = hist2.GetHists()
        sumHist = hist1.Clone("sumHist")
        sumHist.Reset()
        for h in hList:
            sumHist.Add(h)
        retH.Add(sumHist,-1.)
        retH.Divide(sumHist)

        retH.GetXaxis().SetLabelSize(sumHist.GetXaxis().GetLabelSize())
        retH.GetXaxis().SetLabelOffset(sumHist.GetXaxis().GetLabelOffset())
        retH.GetXaxis().SetTitleSize(sumHist.GetXaxis().GetTitleSize())
        retH.GetXaxis().SetTitleOffset(sumHist.GetXaxis().GetTitleOffset())
        retH.GetYaxis().SetLabelSize(sumHist.GetYaxis().GetLabelSize())
        retH.GetYaxis().SetLabelOffset(sumHist.GetYaxis().GetLabelOffset())
        retH.GetYaxis().SetTitleSize(sumHist.GetYaxis().GetTitleSize())
        retH.GetYaxis().SetTitleOffset(sumHist.GetYaxis().GetTitleOffset())

    except(AttributeError):
        #this is the error you get if hist1 is a stack
        print "Did you use a stack as argument 1? please use stack as argument 2!"
        raise AttributeError
    if ymax or ymin:
        retH.GetYaxis().SetRangeUser(0.5,1.5)
        retH.SetLineColor(hist2.GetLineColor())
    ROOT.SetOwnership(retH,0)
    retH.SetYTitle("(data - MC)/MC")
    return retH


def makeDivCan():
    """returns a divided canvas for ratios"""
    Rcanv = ROOT.TCanvas("Rcanv","Rcanv",1024,768)
    Rcanv.cd()
    pad1 = ROOT.TPad("pad1","pad1",0,0.2,1,1)
    pad1.SetNumber(1)
    pad1.SetFillColor(0)
    pad1.SetBorderMode(0)
    pad1.SetFrameFillStyle(0)
    pad1.SetFrameBorderMode(0)
    pad1.SetTickx(0)
    pad1.SetTicky(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.2)
    pad2.SetNumber(2)
    pad2.SetFillColor(0)
    pad2.SetBorderMode(0)
    pad2.SetFrameFillStyle(0)
    pad2.SetFrameBorderMode(0)
    pad2.SetTickx(0)
    pad2.SetTicky(0)


    pad1.Draw()
    pad2.Draw()
    Rcanv.cd()
    #ROOT.SetOwnership(Rcanv,0)
    ROOT.SetOwnership(pad1,0)
    ROOT.SetOwnership(pad2,0)
    #print "check stuff", canv.GetPad(1)
    return Rcanv



def MakeStatProgression(myHisto,title=""):
    """This function returns a function with the statistical precision in each bin"""
    statPrecision = myHisto.Clone()
    ##statPrecision.Reset()
    statPrecision.SetTitle(title)
    for bin in range(myHisto.GetNbinsX()):
        y=statPrecision.GetBinContent(bin);
        err=statPrecision.GetBinError(bin);
        if (y>0): statPrecision.SetBinContent(bin,err/y);
        statPrecision.SetBinError(bin,0);

    return statPrecision
