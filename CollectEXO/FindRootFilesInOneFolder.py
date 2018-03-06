import subprocess
import sys
import fileinput
import os
def getFileList(inDir, inputName):
    tempfile = open(inputName,'w')
    cmd = 'find ' + inDir + ' -mindepth 2 -maxdepth 2'+' -name *_selec*.root'
    #cmd = 'find ' + inDir + ' -mindepth 2 -maxdepth 2'+' -name *Tfile*.root'
    print cmd
    process = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    tempfile.write(process.communicate()[0])
    tempfile.close()
    for line in fileinput.input([inputName], inplace=True):
        sys.stdout.write('root://eoscms/{l}'.format(l=line))
    

inputName_='testscript.txt'
inDir_='/eos/cms/store/group/phys_higgs/HiggsExo/mshi/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/SignalH125a19v1/171025_202130'
#inDir_='/eos/cms/store/group/phys_higgs/HiggsExo/mshi/DoubleMuon/data_RegionE_Thirteen/171024_201451'
#inDir_='/eos/cms/store/group/phys_higgs/HiggsExo/mshi/DoubleMuon/data_RegionF_Thirteen/171024_002045'
getFileList(inDir_,inputName_)

