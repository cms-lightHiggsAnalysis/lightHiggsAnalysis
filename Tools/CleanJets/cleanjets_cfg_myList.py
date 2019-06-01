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
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/r/rhabibul/CMSSW_8_0_17/src/Tools/CleanJets/inputFile_SingleMuon_DYFakeRate.out')

readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
      fileNames = readFiles,
      skipEvents = cms.untracked.uint32(0)
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.CleanJets.outFileName = cms.string('/afs/cern.ch/user/r/rhabibul/CMSSW_8_0_17/src/Tools/CleanJets/OUTPUT.root')

#######################################
# HPS Tau Reconstruction alterations 
#######################################
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.recoTauAK4PFJets08Region.src = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')

process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.ak4PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.combinatoricRecoTaus.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')

#########################################
# Rerunning bTaggin on CleanJets Sample
#########################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")  #Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") #Configuration.StandardSequences.FrontierConditions_CMS.GlobalTag_cff")

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
# Configuring Output
#######################################
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:DIRNAME_NUM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
)

process.p = cms.Path(process.muonsRef*
                     process.CleanJets*
                     process.PFTau*
                          process.pfBTagging)
process.e = cms.EndPath(process.out)
