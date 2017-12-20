import FWCore.ParameterSet.Config as cms
from subprocess import *
import FWCore.Utilities.FileUtils as FileUtils
mylist=FileUtils.loadListFromFile('/afs/cern.ch/work/m/mshi/private/CMSSW_8_0_17/src/CollectEXO/test.txt')
process = cms.Process("SKIM")
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
process.GlobalTag.globaltag = cms.string('80X_dataRun2_2016SeptRepro_v4') # CMSSW 8
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

# b-tagging general configuration
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak4JTA_cff")
# configure the softMuonByIP3d ESProducer and EDProducer
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softLepton_cff import *
from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
process.impactParameterTagInfos=process.impactParameterTagInfos.clone()
process.impactParameterTagInfos.jetTracks = cms.InputTag("ak4JetTracksAssociatorAtVertex")
process.ak4JetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets','ak4PFJetsNoMu','SKIM')
process.ak4JetTracksAssociatorAtVertex.tracks = cms.InputTag("generalTracks")

	
process.btagging = cms.Sequence(
    process.ak4JetTracksAssociatorAtVertex*
    # impact parameters and IP-only algorithms
    process.impactParameterTagInfos*
    (process.trackCountingHighEffBJetTags +
     process.trackCountingHighPurBJetTags +
     process.jetProbabilityBJetTags +
     process.jetBProbabilityBJetTags +
     # SV tag infos depending on IP tag infos, and SV (+IP) based algos
     process.secondaryVertexTagInfos*
     (process.simpleSecondaryVertexHighEffBJetTags +
      process.simpleSecondaryVertexHighPurBJetTags 
#+
     # process.combinedSecondaryVertexBJetTags
     ) +
     process.ghostTrackVertexTagInfos*
     process.ghostTrackBJetTags)##  +
##     process.softPFMuonsTagInfos*
##     process.softPFMuonBJetTags *
##     process.softPFElectronsTagInfos*
##     process.softPFElectronBJetTags
)
#only proceed if event is a true W-->munu event

#require event to fire IsoMu24_eta2p1
process.TriggerAnalyzer0=cms.EDAnalyzer("TriggerAnalyzer")
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     #HLTPaths = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*" ),
     #HLTPaths = cms.vstring("HLT_Mu7p5_Track7_Upsilon_v*", "HLT_Mu7p5_Track3p5_Upsilon_v*","HLT_Mu7p5_Track2_Upsilon_v*", "HLT_Mu7p5_L2Mu2_Upsilon_v*","HLT_Dimuon8_Upsilon_Barrel_v*","HLT_Dimuon13_Upsilon_v*","HLT_Dimuon0_Upsilon_Muon_v*" ),
     HLTPaths = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
)
process.hTozzTo4leptonsMuonCalibrator = cms.EDProducer("HZZ4LeptonsMuonCalibrator",
    muonCollection = cms.InputTag("muons"),
    identifier = cms.string("DATA_80X_13TeV"),
    isData         = cms.bool(True)
)
process.RandomNumberGeneratorService = cms.Service(
    "RandomNumberGeneratorService",
    RochesterCorr = cms.PSet(
        initialSeed = cms.untracked.uint32(3),
        engineName = cms.untracked.string('TRandom3')
    ),
)

process.RochesterCorr=cms.EDProducer("Rochester",
    muonCollection = cms.InputTag("muons"),
    identifier = cms.string("DATA_80X_13TeV"),
    isData         = cms.bool(True),
    initialSeed = cms.untracked.uint32(89),
    engineName = cms.untracked.string('TRandom3'),
    fp=cms.FileInPath("Rochester/RochesterSub/data/rcdata.2016.v3/config.txt")
)   
MU_CUT=("pt>5.0 && abs(eta)<2.4")
process.PreMuons = cms.EDFilter('MuonRefSelector',
                                 src = cms.InputTag("RochesterCorr","RochesterMu","SKIM"),
                                 cut = cms.string(MU_CUT),
                                 filter = cms.bool(True)
)
process.Mu45Selector = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag('MuonIWant'),
    genParticleTag = cms.InputTag('genParticles'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    MatchCut = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu50_v3", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu50_v3","","HLT"),#TTBar background is v2
    #theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q","","HLT2"),#v2
    theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q","","HLT"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1),
    outFileName=cms.string("Mu45Selector.root")
)
process.AllPreMuonsID=cms.EDFilter(
  'MuonsID',
  muonTag=cms.InputTag('PreMuons'),
  vtxTag= cms.InputTag('offlinePrimaryVertices'),
  muonID=cms.string('medium')
)
process.Isolate=cms.EDFilter('CustomDimuonSelector',
                                muonTag=cms.InputTag('AllPreMuonsID'),
                                isoMax=cms.double(0.25),
                                isoMin=cms.double(0.0),
                                baseMuonTag=cms.InputTag('muons'),
                                particleFlow=cms.InputTag('particleFlow'),
                                minNumObjsToPassFilter=cms.uint32(2)
)
process.HighestPtAndMuonSignDRSelector=cms.EDFilter(
                'HighestPtAndMuonSignDRSelector',
                muonTag=cms.InputTag('AllPreMuonsID'),
                dRCut=cms.double(-1),
                passdR=cms.bool(True),
                Mu1PtCut=cms.double(20.0),
                Mu2PtCut=cms.double(20.0),
                oppositeSign = cms.bool(True) # False for SameSignDiMu, True regular
)
process.NonIsolate=cms.EDFilter('CustomDimuonSelector',
				muonTag=cms.InputTag('Mu1Mu2EtaCut'),
				isoMax=cms.double(-1),
				isoMin=cms.double(1.0),
				baseMuonTag=cms.InputTag('muons'),
				particleFlow=cms.InputTag('particleFlow'),
				minNumObjsToPassFilter=cms.uint32(2)
)
process.PtEtaCut = cms.EDFilter('PTETACUT',
                                 muonTag=cms.InputTag('Isolate'),
                                 Eta=cms.double(-1),
                                 Pt=cms.double(0.0),
                                 minNumObjsToPassFilter=cms.uint32(1)
)
                                

process.Mu3=cms.EDFilter('VetoMuon',
  muonTag=cms.InputTag('MuonIWant'),
  vetoMuonTag=cms.InputTag('Mu1Mu2EtaCut'),
  dRCut=cms.double(0.5),
  minNumObjsToPassFilter=cms.uint32(1)
)

process.Mu3ID = cms.EDFilter('CustomMuonSelector',
                                       baseMuonTag = cms.InputTag('muons'),
                                       muonTag = cms.InputTag('Mu3'),
                                       vtxTag = cms.InputTag('offlinePrimaryVertices'),
                                       muonID = cms.string('loose'),
                                       PFIsoMax = cms.double(-1),
                                       detectorIsoMax = cms.double(-1.0),
                                       PUSubtractionCoeff = cms.double(0.5),
                                       usePFIso = cms.bool(True),
                                       passIso = cms.bool(True),
                                       etaMax = cms.double(2.4),
                                       minNumObjsToPassFilter = cms.uint32(1)
                                       )
process.tauMuonPtSelector=cms.EDFilter('PTETACUT',
                                 muonTag=cms.InputTag('Mu3ID'),
                                 Eta=cms.double(2.4),
                                 Pt=cms.double(5.0),
                                 minNumObjsToPassFilter=cms.uint32(1)
)

#clean the jets of soft muons, then rebuild the taus

process.CleanJets.muonSrc=cms.InputTag('tauMuonPtSelector') # 
process.CleanJets.PFCandSrc = cms.InputTag('particleFlow')
process.CleanJets.cutOnGenMatches = cms.bool(False)
process.CleanJets.outFileName = cms.string('CleanJets.root')
process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')
process.recoTauAK4PFJets08Region.src = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')
process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')
process.ak4PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')
process.combinatoricRecoTaus.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')

process.recoTauCommonSequence = cms.Sequence(   process.Mu3*
						process.Mu3ID*
						process.tauMuonPtSelector*
						process.CleanJets*
						process.ak4PFJetTracksAssociatorAtVertex*
						process.recoTauAK4PFJets08Region*
						process.recoTauPileUpVertices*
						process.pfRecoTauTagInfoProducer
)

process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence) # Kyle Changed  This

from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softLepton_cff import *
from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoBTag.Combined.combinedMVA_cff import *
from RecoBTag.CTagging.RecoCTagging_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
process.pfImpactParameterTagInfos.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')
process.softPFMuonsTagInfos.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')
process.softPFElectronsTagInfos.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM')

process.pfBTagging = cms.Sequence(
    (
      # impact parameters and IP-only algorithms
      pfImpactParameterTagInfos *
      ( pfTrackCountingHighEffBJetTags +
        pfJetProbabilityBJetTags +
        pfJetBProbabilityBJetTags +

        # SV tag infos depending on IP tag infos, and SV (+IP) based algos
        pfSecondaryVertexTagInfos *
        ( pfSimpleSecondaryVertexHighEffBJetTags +
          pfCombinedSecondaryVertexV2BJetTags
        )
        + inclusiveCandidateVertexing *
        pfInclusiveSecondaryVertexFinderTagInfos *
        pfSimpleInclusiveSecondaryVertexHighEffBJetTags *
        pfCombinedInclusiveSecondaryVertexV2BJetTags

      ) +

      # soft lepton tag infos and algos
      softPFMuonsTagInfos *
      softPFMuonBJetTags
      + softPFElectronsTagInfos *
      softPFElectronBJetTags
    ) *

    # overall combined taggers
    ( #CSV + soft-lepton + jet probability discriminators combined
      pfCombinedMVAV2BJetTags

    )
)

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets
#this will produce a ref to the cleaned tau collection
process.muHadTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauHadIsoTag = cms.InputTag('hpsPFTauDiscriminationByCombinedIsolationDeltaBetaCorrRaw3Hits', '',
                                'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFindingNewDMs', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets','valMap','SKIM'),
    overlapCandTag = cms.InputTag('Mu45Selector'),
    overlapCandTag1= cms.InputTag('Mu1Mu2EtaCut'),
    passDiscriminator = cms.bool(True),
    pTMin = cms.double(10.0),
    etaMax = cms.double(2.4),
    isoMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1),
    outFileName=cms.string('muHadTauSelector.root')
    )
process.muHadIsoTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauHadIsoTag = cms.InputTag('hpsPFTauDiscriminationByCombinedIsolationDeltaBetaCorrRaw3Hits', '',
                                'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFindingNewDMs', '', 'SKIM'),
    ),
    jetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets','valMap','SKIM'),
    overlapCandTag = cms.InputTag('Mu45Selector','','SKIM'),
    overlapCandTag1=cms.InputTag('Mu1Mu2EtaCut','','SKIM'),
    passDiscriminator = cms.bool(True),
    pTMin=cms.double(10.0),
    etaMax = cms.double(2.4),
    isoMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1),
    outFileName=cms.string('muHadIsoTauSelector.root')
    )

process.muHadNonIsoTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    tauTag = cms.InputTag('muHadTauSelector','','SKIM'),
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauHadIsoTag = cms.InputTag('hpsPFTauDiscriminationByCombinedIsolationDeltaBetaCorrRaw3Hits', '',
                                'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw', '', 'SKIM')
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
process.MassCut=cms.EDFilter('Mu1Mu2MassFilter',
                                    Mu1Mu2=cms.InputTag('HighestPtAndMuonSignDRSelector'),
                                    minMass=cms.double(60.0),
                                    maxMass=cms.double(120.0)
)
process.Mu1Mu2Analyzer=cms.EDAnalyzer(
  'Mu1Mu2Analyzer',
  Mu1Mu2=cms.InputTag("MassCut",'','SKIM'),
  Mu2PtBins=cms.vdouble(x for x in range(0, 200)),
  invMassBins=cms.vdouble(x for x in range(60, 120)),
  MC=cms.bool(False),
  pfMet=cms.InputTag("pfMet","","RECO"),
  fp=cms.FileInPath("GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/data/pileupWeightForEraC.root"),
  PUTag=cms.InputTag("addPileupInfo","","HLT"),
  Generator=cms.InputTag("generator","","SIM")
)
process.GetRunNumber = cms.EDAnalyzer('GetRunNumber')
#sequences
process.MuMuSequenceSelector=cms.Sequence(
        #process.TriggerAnalyzer0*
        process.HLTEle*
        process.RochesterCorr*
        process.PreMuons*
        process.AllPreMuonsID*
        process.HighestPtAndMuonSignDRSelector*
        process.MassCut*
        process.Mu1Mu2Analyzer*
        process.GetRunNumber
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
