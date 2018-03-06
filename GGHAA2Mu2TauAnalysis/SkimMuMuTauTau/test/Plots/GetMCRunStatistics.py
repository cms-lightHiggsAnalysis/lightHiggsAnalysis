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
def getByRunNumber():
    runNumber=[]
    CrossSection=[]
    Count=0;
    with open("2016_eraF.csv") as inf:
        reader=csv.reader(inf, delimiter=',')
        start=False
        for row in reader:
            if len(row)==0:
                continue
            if checkIfStartWithSpecial(row[0]):
                continue
            Count=Count+1
            if Count==1:
                minRow=int(row[0].split(':')[0])
            runNumber.append(int(row[0].split(':')[0])-minRow+1)
            CrossSection.append( float(row[-1]))
        return [runNumber, CrossSection]

def checkIfStartWithSpecial(Str):
    if Str[0]=='#':
        return True
    return False
def getByRunNumberData():
    tfile1=ROOT.TFile.Open("eraF.root")
    data=tfile1.Get("GetRunNumber/RunNumber")
    print(data)
    print(data.GetSize())
    runNumber=[]
    Statistics=[]
    Count=0
    for i in range(data.GetSize()):
        if data.GetBinContent(i)!=0:
            Count=Count+1
            runNumber.append(i)
            Statistics.append(data.GetBinContent(i))
    return [runNumber, Statistics]

MC_ScaleFactor=getByRunNumber()
print("Printing out MC and Data Statistics")
print(MC_ScaleFactor)
Data_Statistics=getByRunNumberData()
print(Data_Statistics)
# read from root file
tfile = ROOT.TFile.Open('DY_Nine.root')
out=ROOT.TFile.Open("HistogramTitle.root", "RECREATE")
 
tree=tfile.Get("lumiTree/LumiTree")
leafValues=array("i", [0])
tree.SetBranchAddress("summedWeights", leafValues)


mc = tfile.Get('Mu1Mu2Analyzer/Eta of Mu2')
test=mc.Integral()
scale1=1.0/(1.925*1e12)*1921.8*3*1.0
f_out=ROOT.TFile("output.root", "RECREATE")
f_out.cd()
c1=ROOT.TCanvas("c1","",800,600)
Intersect=set(MC_ScaleFactor[0]).intersection(Data_Statistics[0])
print(len(Intersect))
runhist=ROOT.TH1D("histo_MC", "histo_MC",len(MC_ScaleFactor[0]),0, len(MC_ScaleFactor[0]))
datahist=ROOT.TH1D("histo_Data", "histo_Data", len(Data_Statistics[0]), 0, len(Data_Statistics[0]))
loopthrough=list(Intersect)
loopthrough.sort()
print(loopthrough)
mc_indx=0
data_indx=0
for i in range(len(loopthrough)):
    while MC_ScaleFactor[0][mc_indx]<loopthrough[i]:
        mc_indx+=1
    while Data_Statistics[0][data_indx]<loopthrough[i]:
        data_indx+=1
    FinalScale=MC_ScaleFactor[1][mc_indx]*test*scale1
    runhist.SetBinContent(i, (FinalScale))
    datahist.SetBinContent(i, Data_Statistics[1][data_indx])
    
f_out.cd()
c1.cd()
runhist.Draw()
f_out.Write()
f_out.Close()

