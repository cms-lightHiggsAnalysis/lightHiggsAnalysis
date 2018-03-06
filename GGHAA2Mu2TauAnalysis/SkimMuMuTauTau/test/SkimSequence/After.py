import FWCore.ParameterSet.Config as cms
from subprocess import *
import FWCore.Utilities.FileUtils as FileUtils
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
mylist=FileUtils.loadListFromFile('/afs/cern.ch/work/m/mshi/private/CMSSW_8_0_17/src/CollectEXO/AllRootFiles/myListAfterEraF.txt')
process = cms.Process("After")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
                SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*mylist))
process.TFileService = cms.Service("TFileService",
    fileName =  cms.string('AfterTFile.root')
)
process.GetRunNumber = cms.EDAnalyzer('GetRunNumber')
process.p = cms.Path(process.GetRunNumber)
