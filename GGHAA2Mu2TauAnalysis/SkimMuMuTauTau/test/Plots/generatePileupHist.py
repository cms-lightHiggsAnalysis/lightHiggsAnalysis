#!/usr/bin/env python
import sys
import os

import ROOT


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    histName = 'pileup'
    fileName = 'pileup.root'
    
    # 80X moriond pileup
    #from SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi import mix
    probValue = [1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05]
    pileupDist = [float(x) for x in probValue]

    rootfile = ROOT.TFile(fileName,'recreate')
    
    # create mc pileup dist
    histmc = ROOT.TH1D(histName+'_MC',histName+'_MC',len(pileupDist),0,len(pileupDist))
    for b,val in enumerate(pileupDist):
        histmc.SetBinContent(b+1,val)
    histmc.Scale(1./histmc.Integral())
    
    histmc.Write()
    
    # read data
    datafile = ROOT.TFile('PileupEraF.root')
    histdata = datafile.Get('pileup')
    histdata.SetTitle("Data")
    histdata.SetName("data")
    histdata.Scale(1./histdata.Integral())
    rootfile.cd()
    histdata.Write()
    
    # now use to get scalefactors
    numbins = min([histdata.GetNbinsX(),histmc.GetNbinsX()])
    histscale = ROOT.TH1D(histName+'_scale',histName+'_scale',numbins,0,numbins)
    for b in range(numbins):
        d = histdata.GetBinContent(b+1)
        m = histmc.GetBinContent(b+1)
        sf = float(d)/m if m else 0.
        histscale.SetBinContent(b+1,sf)
    histscale.Write()
    
    rootfile.Write()
    rootfile.Close()



if __name__ == "__main__":
    status = main()
    sys.exit(status)
                    
