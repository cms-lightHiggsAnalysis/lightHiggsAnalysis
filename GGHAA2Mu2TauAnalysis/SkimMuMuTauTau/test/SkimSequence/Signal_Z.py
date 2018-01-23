import FWCore.ParameterSet.Config as cms
from subprocess import *
import FWCore.Utilities.FileUtils as FileUtils
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
mylist=FileUtils.loadListFromFile('/afs/cern.ch/work/m/mshi/private/CMSSW_8_0_17/src/CollectEXO/DataMini.txt')
process = cms.Process("MINIAODSKIM")
#Debug utils
#process.ProfilerService = cms.Service (
#      "ProfilerService",
#       firstEvent = cms.untracked.int32(2),
#       lastEvent = cms.untracked.int32(500),
#       paths = cms.untracked.vstring('schedule')
#)
#process.SimpleMemoryCheck = cms.Service(
#    "SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
#PDG IDs
A_PDGID = 36
Z_PDGID = 23
W_PDGID = 24
TAU_PDGID = 15
MU_PDGID = 13
NUMU_PDGID = 14
D_PDGID = 1
U_PDGID = 2
S_PDGID = 3
C_PDGID = 4
B_PDGID = 5
T_PDGID = 6
G_PDGID = 21
ANY_PDGID = 0

#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

#tau hadronic decay types
TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15

#no consideration of pT rank
ANY_PT_RANK = -1

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
                SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*mylist))

process.source.inputCommands = cms.untracked.vstring("keep *")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
#for mu-less jets
process.load('Configuration.StandardSequences.MagneticField_cff') #I changed it from: process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') # Kyle Added this
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi') # Kyle Added this
process.GlobalTag.globaltag = cms.string('80X_dataRun2_2016SeptRepro_v7') # CMSSW 8
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")
process.load('Tools/CleanJets/cleanjets_cfi')
process.load('Configuration.StandardSequences.Services_cff')
#define a parameter set to be passed to all modules that utilize GenTauDecayID for signal taus
AMuMuPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

ATauTauPSet=cms.PSet(momPDGID = cms.vint32(A_PDGID),
				   chargedHadronPTMin=cms.double(0.0),
				   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))
# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences. 
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.metTools import *


# to use tau-cleaned jet collection uncomment the following: 
#getattr(process,"pfNoTau"+postfix).enable = True

# to switch default tau to HPS tau uncomment the following: 
#adaptPFTaus(process,"hpsPFTau",postfix=postfix)


#output commands
skimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    "keep *",
    "drop *_*ak7*_*_*",
    "drop *_*ak5*_*_*",
    "drop *_*ak8*_*_*",
    "drop *_*GenJets_*_*",
    "drop *_ak4CaloJets*_*_*",
    "drop *_ak4TrackJets*_*_*",
    "drop *_*jetCentral*_*_*",
    "drop *_hfEMClusters_*_*",
    "drop *_eid*_*_*",
    "drop *_muonMETValueMapProducer_muCorrData_*",
    "drop *_muons_muonShowerInformation_*",
    "drop *_muons_combined_*",
    "drop *_muons_csc_*",
    "drop *_muons_dt_*",
    "drop l1extraL1HFRings_*_*_*",
    "drop *_muonsFromCosmics*_*_*",
    "drop recoCaloClusters_*_*_*",
    "drop recoPreshowerCluster*_*_*_*",
    "drop *_hfRecoEcalCandidate_*_*",
    "drop *_generalV0Candidates_*_*",
    "drop *_selectDigi_*_*",
    "drop *_castorreco_*_*",
    "drop *_reduced*RecHits*_*_*",
    "drop *_PhotonIDProd_*_*",
    "drop *_*_*photons*_*",
    "drop *_dedx*_*_*",
    "drop *_*_cosmicsVeto_*",
    "drop *_muonMETValueMapProducer_*_*",
    "drop *_BeamHaloSummary_*_*",
    "drop *_GlobalHaloData_*_*",
    "drop *_*_uncleanOnly*_*",
    "drop recoCaloMET_*_*_*",
    "drop recoConversion_*_*_*",
    "drop *_CastorTowerReco_*_*",
    "drop *_uncleanedOnlyGsfElectron*_*_*",
    "drop recoJPTJet_*_*_*",
    "drop recoPFMET_*_*_*",
    "drop *_photons_*_*",
    "drop *_photonCore_*_*",
    "drop *_hpsPFTauDiscrimination*_*_RECO",
    "drop *_hpsPFTauProducer_*_RECO",
    "drop *_recoTauAK4PFJets08Region_*_SKIM",
    "drop *_combinatoricRecoTaus_*_SKIM",
    "drop *_hpsPFTauProducerSansRefs_*_SKIM",
    "drop *_pfRecoTauTagInfoProducer_*_SKIM",
    "drop *_recoTauPileUpVertices_*_SKIM",
    "drop *_correctedHybridSuperClusters_*_*",
    "drop *_correctedMulti5x5SuperClustersWithPreshower_*_*",
    "drop *_*phPFIsoValue*04PFIdPFIso_*_*",
    "drop *_*TagInfos*_*_*",
    "drop *_*NoNu_*_*",
    "drop *_clusterSummaryProducer_*_*",
    "drop *_hcalnoise_*_*",
    "drop *_castorDigis_*_*",
    "drop *_hcalDigis_*_*",
    "drop *_tevMuons_*_*",
    "drop *_logErrorHarvester_*_*",
    "drop *_particleFlowTmp_*_*",
    "drop *_particleFlowRecHit*_*_*"
    )
  ) 

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(process,
#                           isData=True, # (or False),
#                           jetCollUnskimmed = "updatedPatJetsUpdatedJEC"
#                           )
from RecoMET.METFilters.metFilters_cff import HBHENoiseFilterResultProducer, HBHENoiseFilter, HBHENoiseIsoFilter, hcalLaserEventFilter
from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter, eeBadScFilter, ecalLaserCorrFilter, EcalDeadCellBoundaryEnergyFilter
from RecoMET.METFilters.metFilters_cff import primaryVertexFilter, CSCTightHaloFilter, CSCTightHaloTrkMuUnvetoFilter, CSCTightHalo2015Filter, globalTightHalo2016Filter, globalSuperTightHalo2016Filter, HcalStripHaloFilter
from RecoMET.METFilters.metFilters_cff import trackingFailureFilter, trkPOGFilters, manystripclus53X, toomanystripclus53X, logErrorTooManyClusters
from RecoMET.METFilters.metFilters_cff import chargedHadronTrackResolutionFilter, muonBadTrackFilter
from RecoMET.METFilters.metFilters_cff import BadChargedCandidateSummer16Filter, BadPFMuonSummer16Filter #2016 ICHEP version
from RecoMET.METFilters.metFilters_cff import metFilters

Flag_goodVertices= cms.Path(primaryVertexFilter)
Flag_globalTightHalo2016Filter = cms.Path(globalTightHalo2016Filter)
Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilterResultProducer * HBHENoiseFilter)
Flag_HBHENoiseIsoFilter = cms.Path(HBHENoiseFilterResultProducer * HBHENoiseIsoFilter)
Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)

process.METTrigger =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","PAT"),
     HLTPaths = cms.vstring("Flag_goodVertices","Flag_globalTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","EcalDeadCellTriggerPrimitiveFilter"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(False), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
)
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#require event to fire IsoMu24_eta2p1
process.TriggerAnalyzer0=cms.EDAnalyzer("TriggerAnalyzer")
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
)
process.RandomNumberGeneratorService = cms.Service(
    "RandomNumberGeneratorService",
    RochesterCorr = cms.PSet(
        initialSeed = cms.untracked.uint32(3),
        engineName = cms.untracked.string('TRandom3')
    ),
)

process.RochesterCorr=cms.EDProducer("Rochester",
    muonCollection = cms.InputTag("slimmedMuons"),
    identifier = cms.string("DATA_80X_13TeV"),
    isData         = cms.bool(True),
    initialSeed = cms.untracked.uint32(89),
    engineName = cms.untracked.string('TRandom3'),
    fp=cms.FileInPath("Rochester/RochesterSub/data/rcdata.2016.v3/config.txt")
)   
process.PreMuons = cms.EDFilter('PTETACUT',
            muonTag=cms.InputTag("RochesterCorr","RochesterMu"),
            Eta=cms.double(2.4),
            Pt=cms.double(5.0),
            minNumObjsToPassFilter=cms.uint32(2)
)
process.MuonsIDdxydz=cms.EDFilter(
  'MuonsID',
  muonTag = cms.InputTag('PreMuons'),
  muonID = cms.string('medium')
)

process.TrigMuMatcher=cms.EDFilter(
        'TrigMuMatcher',
        muonsTag = cms.InputTag('MuonsIDdxydz'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        prescales = cms.InputTag("patTrigger"),
        triggerObjects = cms.InputTag("selectedPatTrigger"),
        trigNames = cms.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v"),
        dRCut = cms.double(.15),
        mu1PtCut = cms.double(26.0)
)
process.GetMuOne = cms.EDFilter(
       'GetHighestPt',
       muonTag = cms.InputTag('TrigMuMatcher')
)
process.Mu2Iso=cms.EDFilter(
  'Mu2Iso',
  muonTag = cms.InputTag('MuonsIDdxydz'),
  mu1Tag = cms.InputTag('GetMuOne'),
  relIsoCutVal = cms.double(0.25), # .25 for iso, -1 for ignoring iso
  passRelIso = cms.bool(True) #False = Non-Iso DiMu, True = Iso-DiMu
)

process.DiMuSigndRSelector=cms.EDFilter(
                'DiMuSigndRSelector',
                mu1Tag = cms.InputTag('GetMuOne'),
                muonsTag = cms.InputTag('Mu2Iso'),
                dRCut = cms.double(-1),
                passdR = cms.bool(True),
                oppositeSign = cms.bool(True) # False for SameSignDiMu, True regular
)

process.GetMuTwo = cms.EDFilter(
       'GetHighestPt',
       muonTag = cms.InputTag('DiMuSigndRSelector')
)

process.Mu1Mu2 = cms.EDFilter(
        'CombineMu1Mu2',
        mu1Tag = cms.InputTag('GetMuOne'),
        mu2Tag = cms.InputTag('GetMuTwo')
)
process.MassCut=cms.EDFilter('Mu1Mu2MassFilter',
                                    Mu1Mu2=cms.InputTag('Mu1Mu2'),
                                    minMass=cms.double(60.0),
                                    maxMass=cms.double(120.0)
)
process.Mu1Mu2Analyzer=cms.EDAnalyzer(
  'Mu1Mu2Analyzer',
  Mu1Mu2=cms.InputTag("MassCut",'','MINIAODSKIM'),
  Mu1=cms.InputTag("GetMuOne"),
  Mu2PtBins=cms.vdouble(x for x in range(0, 200)),
  invMassBins=cms.vdouble(x for x in range(60, 120)),
  MC=cms.bool(False),
  pfMet=cms.InputTag("slimmedMETs"),
  fp=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/pileup.root"),
  fpIDs_BToF=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/EfficienciesAndSF_BCDEF.root"),
  fpIDs_GH=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/EfficienciesAndSF_GH.root"),
  fpISOs_BToF=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/EfficienciesAndSF_BCDEF_ISO.root"),
  fpISOs_GH=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/EfficienciesAndSF_GH_ISO.root"),
  fpTrack=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/Tracking_EfficienciesAndSF_BCDEFGH.root"),
  fpTrigger_BToF=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/EfficienciesAndSF_RunBtoF.root"),
  fpTrigger_GH=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/EfficienciesAndSF_Period4.root"),
  PUTag=cms.InputTag("addPileupInfo","","HLT"),
  Generator=cms.InputTag("generator"),
  BadChargedCandidateFilter=cms.InputTag("BadChargedCandidateFilter"),
  BadPFMuonFilter=cms.InputTag("BadPFMuonFilter")
)
process.GetRunNumber = cms.EDAnalyzer('GetRunNumber')
#sequences
process.MuMuSequenceSelector=cms.Sequence(
        #process.TriggerAnalyzer0*
        process.jecSequence*
        process.METTrigger*
        process.HLTEle*
        process.RochesterCorr*
        process.PreMuons*
        process.MuonsIDdxydz*
        process.TrigMuMatcher*
        process.GetMuOne*
        process.Mu2Iso*
        process.DiMuSigndRSelector*
        process.GetMuTwo*
        process.Mu1Mu2*
        process.MassCut*
        process.BadPFMuonFilter *
        process.BadChargedCandidateFilter *
        process.Mu1Mu2Analyzer
)


process.antiSelectionSequence = cms.Sequence(process.MuMuSequenceSelector
)


process.antiSelectedOutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = skimEventContent.outputCommands,
    fileName = cms.untracked.string('RegionB_selection.root')
    )
#sequences
process.TFileService = cms.Service("TFileService",
    fileName =  cms.string('TFiledataZ.root')
)
#no selection path
process.p = cms.Path(process.antiSelectionSequence)
process.e = cms.EndPath(process.antiSelectedOutput)
