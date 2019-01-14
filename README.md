# lightHiggsAnalysis

This analysis repository is a search for a light higgs particle could potentially be found using Large Hadron Collider data collected by CERN. The following instruction is a guide to install this package on your working computer correctly. I will introduce function of specific packages after the instruction.

Go to https://github.com/cms-lightHiggsAnalysis/lightHiggsAnalysis.git and click "Fork" in the upper right-hand corner

Execute:

git clone https://github.com/cms-lightHiggsAnalysis/lightHiggsAnalysis.git

cp -r lightHiggsAnalysis/* .

scram b -j 16

***********************************************************************

## Folders are organized as following: 

### GGHAA2Mu2TauAnalysis 

Contains our major analysis code 

###SkimSequence folder

$CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/ is where our selection sequence developed. It can skim officially generated MC background down to order 10^-6 while keep most of signal.    


And pleas notice, all the variables can be defined by user. 

di-mu isolation

process.Isolate=cms.EDFilter('CustomDimuonSelector',

                                muonTag=cms.InputTag('Mu1Mu2EtaCut'),
                                
                                isoMax=cms.double(1.0),
                                
                                isoMin=cms.double(0.0),
                                
                                baseMuonTag=cms.InputTag('muons'),
                                
                                particleFlow=cms.InputTag('particleFlow'),
                                
                                minNumObjsToPassFilter=cms.uint32(2)
)

process.muHadNonIsoTauSelector = cms.EDFilter(

    'CustomTauSepFromMuonSelector',
    
    tauTag = cms.InputTag('muHadTauSelector','','SKIM'),
    
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    
    tauHadIsoTag = cms.InputTag('hpsPFTauDiscriminationByCombinedIsolationDeltaBetaCorrRaw3Hits', '',
    
                                'SKIM'),
                                
    tauDiscriminatorTags = cms.VInputTag(
    
    cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT', '', 'SKIM')
    
    ),
    
    jetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM'),
    
    muonRemovalDecisionTag = cms.InputTag('CleanJets','valMap','SKIM'),
    
    overlapCandTag = cms.InputTag('Mu45Selector'),
    
    overlapCandTag1=cms.InputTag('Mu1Mu2EtaCut','','SKIM'),
    
    passDiscriminator = cms.bool(False),
    
    pTMin=cms.double(10.0),
    
    etaMax = cms.double(2.4),
    
    isoMax = cms.double(-1.0),
    
    dR = cms.double(0.5),
    
    minNumObjsToPassFilter = cms.uint32(1),
    
    outFileName=cms.string('muHadNoIsoTauSelector.root')
    
)

Specifically, change passDiscriminator from "False" to "True" will yield isolatioon instead of anti-isolation. You can also speficy isoMax, and default is -1, which is no limit.



Detailed information to get this piece of selection work: type

cd $CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/

cmsRun RegionBSkim.py

You will obtain a .root file that stores events after this selection sequence. In particular, this RegionBSkim.py gives you selection sequence of getting B region data defined as "isolated di-mu and non-isolated di-tau selection plus passing general selection sequence".

### AMuTriggerAnalyzer and MuMuTauTauRecoAnalyzer 

Contains analyzers targeting for analyze di-muon behavior.

### Plotting

we have plotting scripts that can give us correctly normalized plots of stack of all background. How to use it is that you cd into $CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/TauAnalyzer/test

and type into following command

 ./runStandardAnalysis.sh 1

It is still under development. It aims to produce more maybe all plots that we are interested in just one go.  

### QCDBackground
Specifically analyze MC QCD background


### CollectEXO

It contains all the files of background and signal Monte Carlo after selection sequence. These MC results are stored on eos space of higgs Exotic group.

### GetStatistics

Everytime you submit job to crab, it divide one single job to smaller jobs, and since selection sequence table is written in log files, this scirpt collects all these tables from your entire job area, and calculate the final statistics. 

You need to change https://github.com/MengyaoShi/lightHiggsAnalysis/blob/master/GetStatistics/Files.txt to sepecify where your log files are store at. Currently it focus on getting selection table for QCD Monte Carlo. But it's flexible and you can study anything you are interested in.

### Tools 

package has a separate readme https://github.com/MengyaoShi/lightHiggsAnalysis/blob/master/Tools/README.md It contains packages has more general functions support this analysis. 


## Workflow

a) In $CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence, define your selection sequence

b) In $CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/crabConfig, define corresponding crab submission configuration file, and submit job to crab server.

c) cd $CMSSW_BASE/src/CollectEXO/ 

mkdir AllRootFiles

Use $CMSSW_BASE/src/CollectEXO/FindRootFiles.py to get collection of all .txt files that contains root file path of skimmed results.

It will look like, in CollectEXO/AllRootFiles

DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt

DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt

QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8.txt

SingleMuon.txt

SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8.txt

SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8.txt

SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-750_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8.txt

SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8.txt

TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt

For example, in TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt,

It contains,

root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/mshi/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTBar_RegionBWithMassCut/170416_022149/0000/RegionB_selection_1.root

root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/mshi/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTBar_RegionBWithMassCut/170416_022149/0000/RegionB_selection_10.root

...

170416_022149 is the crab submission job number.

d) In $CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/QCDBackground, there is a Summary.sh

./Summary.sh, you will submit analysis of each background/signal in lxplus. Since we want plot on same physics measurement, for example pt of highest pt muon, missing energy, and so on, we run similar code on different backgrounds/signals, I wrote a template named test_skimmed.py that be fed into Summary.sh, and copy this python file 10 times and change a bit according to which backgrounds we analyze. In the end of Summary.sh, it submit jobs with all modified versions of analysis template using 'bsub', and start to analyze all backgrounds and signals. Wait until you get .root files, use "bjobs" to check status.

.root files will contains histograms. Modify test_skimmed.py to include new plots.

e) In StandardPlotFormat I use several simple scripts to stack all backgrounds together and compare with signal. Notice this is a temporary solution. In $CMSSW_BASE/src/GGHAA2Mu2TauAnalysis/Plotting/test as mentioned in above sections, it's a non-temporary solution that's under development. 
