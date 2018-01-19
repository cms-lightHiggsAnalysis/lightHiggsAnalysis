import os
import sys
import math

import ROOT

import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
from array import array
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPalette(1)
#GEt byRun number
import csv
# arguments
savename = 'zmass'
lumi = 34872.0
xaxis = 'Pt of Mu1(GeV)'
yaxis = 'Events Statistics(EraBCDEFGH 2016)'
logy = True
logx = False
ymin = 1
ymax = None
numcol = 1
legendpos = 33
lumipos = 11
isprelim = True
ratiomin = 0.8
ratiomax = 1.2
rangex = [0,300]
Luminosity=34872.0
# read from root file
tfile = ROOT.TFile.Open('./RootFiles/AllData.root')

out=ROOT.TFile.Open("HistogramTitle.root", "RECREATE")

#tree=tfile.Get("lumiTree/LumiTree")
#leafValues=array("i", [0])
#tree.SetBranchAddress("summedWeights", leafValues)

tdirectory=tfile.Get("Mu1Mu2Analyzer")
tlist_keys=tdirectory.GetListOfKeys()
print(tdirectory)

hash_={"pt of Mu1 Mu2 (H750a09)":["pt of reco muon", [0,300]], 
"dRMetMu1": ["dRMetMu1",[0,5.0]], 
"MetPt":["MetPt",[0,300]], 
"invMass of Mu1 Mu2 (H750a09)":["invMass of Mu1 Mu2", [0,300]], 
"Mu1Mu2Pt":["Mu1Mu2Pt",[0,300]],
"Mu1Mu2Eta":["Mu1Mu2Eta",[-2.5, 2.5]],
"Mu1PtMu2Pt":["Mu1PtMu2Pt",[0,300]],
"Pt of RecoMu2": ["pt of Mu2",[0,300]],
"pt of Mu1": ["pt of Mu1",[0,300]],
"dRMu1Mu2":["dRMu1Mu2",[0,5.0]],
"dRMu1Mu2Wider":["dRMu1Mu2Wider", [0,5.0]],
"dRVsMu2Pt":["dRVsMu2Pt",[0,4.0]],
"Eta of Mu1": ["Eta of Mu1",[-2.5, 2.5]],
"Eta of Mu2": ["Eta of Mu2",[-2.5, 2.5]],
"NumVertices": ["NumVertices",[0,20]]
}
interested=["Pt of RecoMu2", "invMass of Mu1 Mu2 (H750a09)", "dRMetMu1","MetPt","Mu1Mu2Pt","Mu1Mu2Eta","Eta of Mu1","Eta of Mu2", "pt of Mu1","Mu1PtMu2Pt", "dRMu1Mu2Wider"]

def draw_one_hist(data,  stack, xaxis, rangex):
	data.SetTitle('Observed')
	data.SetMarkerSize(1)
	data.SetMarkerStyle(20)
	# create canvas
	canvas = ROOT.TCanvas(savename,savename,50,50,800,800)
	plotpad = ROOT.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
	plotpad.SetBottomMargin(0.04)
	plotpad.SetRightMargin(0.03)
	out.cd()
	plotpad.Draw()
	plotpad.SetLogy(logy)
	plotpad.SetLogx(logx)
	ratiopad = ROOT.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
	ratiopad.SetTopMargin(0.06)
	ratiopad.SetRightMargin(0.03)
	ratiopad.SetBottomMargin(0.5)
	ratiopad.SetLeftMargin(0.16)
	ratiopad.SetTickx(1)
	ratiopad.SetTicky(1)
	ratiopad.Draw()
	ratiopad.SetLogx(logx)
	plotpad.cd()

	# now draw them
	stack.Draw("hist")
	stack.GetXaxis().SetTitle(xaxis)
	stack.GetYaxis().SetTitle(yaxis)
	stack.GetYaxis().SetLabelSize(0.05)
	if ymax!=None: stack.SetMaximum(ymax)
	if ymin!=None: stack.SetMinimum(ymin)
	if rangex: stack.GetXaxis().SetRangeUser(*rangex)

	# staterr for stack
	stack.GetHistogram().GetXaxis().SetLabelOffset(999)
	staterr = stack.GetStack().Last().Clone("{0}_staterr".format(stack.GetName))
	staterr.SetFillColor(ROOT.kGray+3)
	staterr.SetLineColor(ROOT.kGray+3)
	staterr.SetLineWidth(0)
	staterr.SetMarkerSize(1)
	staterr.SetFillStyle(3013)
	staterr.Draw('e3 same')
	data.Draw('ex0 same')#ex0

	# get the legend
	entries = []
	for hist in reversed(stack.GetHists()):
	    entries += [[hist,hist.GetTitle(),'f']]
	entries += [[data,data.GetTitle(),'ep']]
	#legend = getLegend(*entries,numcol=numcol,position=legendpos)
	legend = ROOT.TLegend(0.7,0.7,0.95,0.95,'','NDC')
	if numcol>1: legend.SetNColumns(int(numcol))
	legend.SetTextFont(42)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	for entry in entries:
	    legend.AddEntry(*entry)
	legend.Draw()

	# cms lumi styling
	period_int = 4
	if plotpad != ROOT.TVirtualPad.Pad(): plotpad.cd()
	CMS_lumi.writeExtraText = isprelim
	CMS_lumi.extraText = "Preliminary"
	CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (float(lumi)/1000.)
	if lumi < 1000:
	    CMS_lumi.lumi_13TeV = "%0.1f pb^{-1}" % (float(lumi))
	CMS_lumi.CMS_lumi(plotpad,period_int,lumipos)


	# the ratio portion for stack
	denom = stack.GetStack().Last().Clone('denom')
	ratiostaterr = denom.Clone("{0}_ratiostaterr".format(denom.GetName))
	ratiostaterr.SetStats(0)
	ratiostaterr.SetTitle("")
	ratiostaterr.GetYaxis().SetTitle("Obs / MC")
	ratiostaterr.SetMaximum(ratiomax)
	ratiostaterr.SetMinimum(ratiomin)
	ratiostaterr.SetMarkerSize(0)
	ratiostaterr.SetFillColor(ROOT.kGray+3)
	ratiostaterr.SetFillStyle(3013)
	ratiostaterr.GetXaxis().SetLabelSize(0.19)
	ratiostaterr.GetXaxis().SetTitleSize(0.21)
	ratiostaterr.GetXaxis().SetTitleOffset(1.0)
	ratiostaterr.GetXaxis().SetLabelOffset(0.03)
	ratiostaterr.GetYaxis().SetLabelSize(0.19)
	ratiostaterr.GetYaxis().SetLabelOffset(0.006)
	ratiostaterr.GetYaxis().SetTitleSize(0.21)
	ratiostaterr.GetYaxis().SetTitleOffset(0.35)
	ratiostaterr.GetYaxis().SetNdivisions(503)
	# bin by bin errors
	for i in range(denom.GetNbinsX()+2):
	    ratiostaterr.SetBinContent(i, 1.0)
	    if denom.GetBinContent(i)>1e-6:  # not empty
		binerror = denom.GetBinError(i) / denom.GetBinContent(i)
		ratiostaterr.SetBinError(i, binerror)
	    else:
		ratiostaterr.SetBinError(i, 999.)
	ratiostaterr.SetXTitle(xaxis)

	unityargs = [rangex[0],1,rangex[1],1] if rangex else [denom.GetXaxis().GetXmin(),1,denom.GetXaxis().GetXmax(),1]
        
	ratiounity = ROOT.TLine(*unityargs)
	ratiounity.SetLineStyle(2)

	# ratio for data
	dataratio = data.Clone('ratio_{0}'.format(data.GetName()))
	for b in range(data.GetNbinsX()):
	    nVal = data.GetBinContent(b+1)
	    nErr = data.GetBinError(b+1)
	    dVal = denom.GetBinContent(b+1)
	    if dVal>1e-6:
		val = nVal/dVal
		err = nErr/dVal
	    else:
		val = 0
		err = 0
	    dataratio.SetBinContent(b+1,val)
	    dataratio.SetBinError(b+1,err)

	# and draw ratio
	if ratiopad != ROOT.TVirtualPad.Pad(): ratiopad.cd()
	ratiostaterr.Draw("e2")
	if rangex: ratiostaterr.GetXaxis().SetRangeUser(*rangex)
	ratiounity.Draw('same')
	dataratio.Draw('ex0 same')

	canvas.Print(xaxis+".png")
        canvas.Close()
for key in tlist_keys:
    if key.GetTitle() in interested:
      print(hash_[key.GetTitle()][0])
      data = tfile.Get('Mu1Mu2Analyzer/'+hash_[key.GetTitle()][0])
      txtfile=open("InputZplot.txt", "r")
      stack=ROOT.THStack('stack', 'stack')
      files = []
      for line in txtfile:
        lis=line.rstrip().split(',')
        print lis
        crossSection=float(lis[1])
        files.append( ROOT.TFile.Open('./RootFiles/'+lis[0]))
        mc=files[-1].Get("Mu1Mu2Analyzer/"+hash_[key.GetTitle()][0])
        MCtree=files[-1].Get("lumiTree/LumiTree")
        summedWeights=0.0
        for entry in MCtree:
	    summedWeights+=MCtree.summedWeights
        print summedWeights
        mc.Scale(1.0/summedWeights*crossSection*Luminosity)
        mc.SetTitle(lis[2])
        mc.SetFillColor(eval(lis[3]))
        stack.Add(mc)
      xaxis=hash_[key.GetTitle()][0]
      myrange=hash_[key.GetTitle()][1]
      draw_one_hist(data, stack,xaxis, myrange)
out.Close()
