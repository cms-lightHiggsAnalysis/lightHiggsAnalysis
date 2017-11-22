import FWCore.ParameterSet.Config as cms

process = cms.Process("CLEANJETS")

#######################################
# Loading all processes and functions
#######################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_v14')
process.load("Tools.CleanJets.cleanjets_cfi")


#######################################
# Declaring Input and configurations
#######################################
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/SkimDir/CMSSW_8_0_17/src/Tools/CleanJets/LessDroppings_NoIsoDiTau_JUL19_B_v3_MostUndropped.out')

readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.CleanJets.outFileName = cms.string('/afs/cern.ch/user/k/ktos/SkimDir/CMSSW_8_0_17/src/Tools/CleanJets/BSUB/DIRNAME/OUTPUT.root')

#########################################
# Rerunning bTaggin on CleanJets Sample
#########################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")  #Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softLepton_cff import *
from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoBTag.Combined.combinedMVA_cff import *
from RecoBTag.CTagging.RecoCTagging_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
process.pfImpactParameterTagInfos.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.softPFMuonsTagInfos.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.softPFElectronsTagInfos.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')

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

#######################################
# HPS Tau Reconstruction alterations 
#######################################
process.Mu3=cms.EDFilter('VetoMuon',
  muonTag=cms.InputTag('MuonIWant'),
  vetoMuonTag=cms.InputTag('Mu1Mu2EtaCut'),
  dRCut=cms.double(0.5),
  genParticleTag = cms.InputTag('genParticles'),
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

process.load('Configuration.StandardSequences.Services_cff')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadPFMuonSummer16Filter_cfi')

## for miniAOD running
process.BadPFMuonFilter.muons = cms.InputTag("muons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("particleFlow")

## for miniAOD running ICHEP
process.BadPFMuonSummer16Filter.muons = cms.InputTag("muons")
process.BadPFMuonSummer16Filter.PFCandidates = cms.InputTag("particleFlow")

process.BadPFMuonFilter.debug = cms.bool(True)
process.BadPFMuonSummer16Filter.debug = cms.bool(True)

#########################
process.load('Configuration.StandardSequences.Services_cff')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateSummer16Filter_cfi')

## for miniAOD running
process.BadChargedCandidateFilter.muons = cms.InputTag("muons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("particleFlow")

process.BadChargedCandidateFilter.debug = cms.bool(True)


process.recoTauCommonSequence = cms.Sequence(   process.Mu3*
                                                process.Mu3ID*
                                                process.tauMuonPtSelector*
                                                process.BadPFMuonFilter*
                                                process.BadChargedCandidateFilter*
                                                process.CleanJets*
                                                process.ak4PFJetTracksAssociatorAtVertex*
                                                process.recoTauAK4PFJets08Region*
                                                process.recoTauPileUpVertices*
                                                process.pfRecoTauTagInfoProducer
)

process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence) # Kyle Changed  This

#######################################
# Configuring Output
#######################################
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/afs/cern.ch/user/k/ktos/SkimDir/CMSSW_8_0_17/src/Tools/CleanJets/BSUB/DIRNAME/DIRNAME_NUM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
)

process.p = cms.Path(process.muonsRef*
                     process.PFTau#*
	   	     #process.pfBTagging
)
process.e = cms.EndPath(process.out)

